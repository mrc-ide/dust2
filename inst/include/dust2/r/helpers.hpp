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

template <typename Iter>
cpp11::sexp export_array_n(Iter iter, std::initializer_list<size_t> dims) {
  const auto len =
    std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<>{});
  cpp11::writable::doubles ret(len);
  std::copy_n(iter, len, ret.begin());
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
