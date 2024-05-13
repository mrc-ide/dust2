#pragma once

#include <numeric>
#include <vector>
#include <dust2/common.hpp>
#include <dust2/cpu.hpp>

namespace dust2 {
namespace r {

inline void check_scalar(cpp11::sexp x, const char * name) {
  if (LENGTH(x) != 1) {
    cpp11::stop("'%s' must be a scalar", name);
  }
}

inline void check_length(cpp11::sexp x, int len, const char * name) {
  if (LENGTH(x) != len) {
    cpp11::stop("'%s' must have length %d", name, len);
  }
}

inline double to_double(cpp11::sexp x, const char * name) {
  check_scalar(x, name);
  if (TYPEOF(x) == REALSXP) {
    return REAL(x)[0];
  }
  if (TYPEOF(x) == INTSXP) {
    return INTEGER(x)[0];
  }
  cpp11::stop("'%s' must be scalar numeric", name);
}

template <typename real_type>
inline std::vector<real_type> to_vector_real(cpp11::sexp x, const char * name) {
  if (TYPEOF(x) == REALSXP) {
    auto x_dbl = cpp11::as_cpp<cpp11::doubles>(x);
    std::vector<real_type> ret(x_dbl.size());
    std::copy(x_dbl.begin(), x_dbl.end(), ret.begin());
    return ret;
  }
  if (TYPEOF(x) == INTSXP) {
    auto x_int = cpp11::as_cpp<cpp11::integers>(x);
    std::vector<real_type> ret(x_int.size());
    std::copy(x_int.begin(), x_int.end(), ret.begin());
    return ret;
  }
  cpp11::stop("'%s' must be a numeric vector", name);
}

inline int to_int(cpp11::sexp x, const char * name) {
  check_scalar(x, name);
  if (TYPEOF(x) == INTSXP) {
    return INTEGER(x)[0];
  }

  if (TYPEOF(x) == REALSXP) {
    double xv = REAL(x)[0];
    if (!cpp11::is_convertible_without_loss_to_integer(xv)) {
      cpp11::stop("'%s' must be integer-like", name);
    }
    return xv;
  }
  cpp11::stop("'%s' must be scalar integer", name);
}

inline size_t to_size(cpp11::sexp x, const char * name) {
  int ret = to_int(x, name);
  if (ret < 0) {
    cpp11::stop("'%s' must be non-negative", name);
  }
  return ret;
}

inline bool to_bool(cpp11::sexp x, const char * name) {
  check_scalar(x, name);
  if (TYPEOF(x) == LGLSXP) {
    return INTEGER(x)[0];
  }
  cpp11::stop("'%s' must be scalar logical", name);
}

template <typename T>
bool is_integer_like(T x, T eps) {
  return std::abs(x - std::round(x)) <= eps;
}

inline double check_time(cpp11::sexp r_time, const char * name) {
  const auto time = to_double(r_time, name);
  const auto eps = 1e-8;
  // We can relax this later and carefully align time onto a grid
  if (!is_integer_like(time, eps)) {
    cpp11::stop("Expected '%s' to be integer-like", name);
  }
  return time;
}

inline double check_dt(cpp11::sexp r_dt) {
  const auto dt = to_double(r_dt, "dt");
  const auto eps = 1e-8;
  if (dt <= 0) {
    cpp11::stop("Expected 'dt' to be greater than 0");
  }
  if (dt > 1) {
    cpp11::stop("Expected 'dt' to be at most 1");
  }
  const auto inv_dt = 1 / dt;
  if (!is_integer_like(inv_dt, eps)) {
    cpp11::stop("Expected 'dt' to be the inverse of an integer");
  }
  return dt;
}

template <typename real_type>
std::vector<real_type> check_time_sequence(real_type time_start,
                                           cpp11::sexp r_time,
                                           const char * name) {
  auto time = to_vector_real<real_type>(r_time, name);
  auto prev = time_start;
  const auto eps = 1e-8;
  for (size_t i = 0; i < time.size(); ++i) {
    const auto t = time[i];
    if (!is_integer_like(t, eps)) {
      cpp11::stop("Expected 'time[%d]' to be integer-like", i + 1);
    }
    if (t <= prev) {
      cpp11::stop("Expected 'time[%d]' (%d) to be larger than the previous value (%d)",
                  i + 1, static_cast<int>(prev), static_cast<int>(t));
    }
    prev = t;
  }
  return time;
}

// The initializer_list is a type-safe variadic-like approach.
inline void set_array_dims(cpp11::sexp data,
                           std::initializer_list<size_t> dims) {
  cpp11::writable::integers r_dim(dims.size());
  auto dim_i = dims.begin();
  for (size_t i = 0; i < dims.size(); ++i, ++dim_i) {
    r_dim[i] = *dim_i;
  }
  data.attr("dim") = r_dim;
}

template <typename It>
cpp11::sexp export_array_n(It it, std::initializer_list<size_t> dims) {
  const auto len =
    std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<>{});
  cpp11::writable::doubles ret(len);
  std::copy_n(it, len, ret.begin());
  set_array_dims(ret, dims);
  return ret;
}

template <typename T>
std::vector<typename T::shared_state> build_shared(cpp11::list r_pars,
                                                   size_t n_groups) {
  const auto grouped = n_groups > 0;
  std::vector<typename T::shared_state> shared;
  if (grouped) {

    if (r_pars.size() != static_cast<int>(n_groups)) {
      cpp11::stop("Expected 'pars' to have length %d to match 'n_groups'",
                  static_cast<int>(n_groups));
    }
    size_t size = 0;
    for (size_t i = 0; i < n_groups; ++i) {
      shared.push_back(T::build_shared(r_pars[i]));
      const auto size_i = T::size(shared[i]);
      if (i == 0) {
        size = size_i;
      } else if (size_i != size) {
        cpp11::stop("Expected state length for group %d to be %d, but it was %d",
                    i + 1, size, size_i);
      }
    }
  } else {
    shared.push_back(T::build_shared(r_pars));
  }
  return shared;
}

template <typename T>
std::vector<typename T::internal_state> build_internal(std::vector<typename T::shared_state> shared) {
  std::vector<typename T::internal_state> internal;
  internal.reserve(shared.size());
  for (auto& s : shared) {
    internal.push_back(T::build_internal(s));
  }
  return internal;
}

template <typename T>
std::vector<typename T::data_type> check_data(cpp11::list r_data,
                                              size_t n_time,
                                              size_t n_groups,
                                              const char * name) {
  const bool grouped = n_groups > 0;
  std::vector<typename T::data_type> data;

  check_length(r_data, n_time, name);

  if (grouped) {
    // There are two ways we might recieve things; as a list-of-lists
    // or as a list matrix.  We might also want to cope with a
    // data.frame but we can probably do that on the R side, and might
    // provide helpers there that throw much nicer errors than we can
    // throw here, really.
    for (size_t i = 0; i < n_time; ++i) {
      auto r_data_i = cpp11::as_cpp<cpp11::list>(r_data[i]);
      check_length(r_data_i, n_groups, "data[i]"); // can do better with sstream
      for (size_t j = 0; j < n_groups; ++j) {
        data.push_back(T::build_data(r_data_i[j]));
      }
    }
  } else {
    for (size_t i = 0; i < n_time; ++i) {
      data.push_back(T::build_data(r_data[i]));
    }
  }

  return data;
}

template <typename T>
void update_pars(dust_cpu<T>& obj, cpp11::list r_pars, bool grouped) {
  if (grouped) {
    const auto n_groups = obj.n_groups();
    if (r_pars.size() != static_cast<int>(n_groups)) {
      cpp11::stop("Expected 'pars' to have length %d to match 'n_groups'",
                  static_cast<int>(n_groups));
    }
    for (size_t i = 0; i < n_groups; ++i) {
      cpp11::list r_pars_i = r_pars[i];
      obj.update_shared(i, [&] (auto& shared) {
                             T::update_shared(r_pars_i, shared);
                            });
    }
  } else {
    obj.update_shared(0, [&] (auto& shared) {
                           T::update_shared(r_pars, shared);
                          });
  }
}

}
}
