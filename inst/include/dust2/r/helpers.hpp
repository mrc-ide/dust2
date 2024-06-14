#pragma once

#include <cstring>
#include <numeric>
#include <vector>
#include <dust2/common.hpp>
#include <cpp11.hpp>

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

inline cpp11::integers as_integers(cpp11::doubles x, const char * name) {
  const auto len = x.size();
  cpp11::writable::integers ret(x.size());
  for (auto i = 0; i < len; ++i) {
    if (!cpp11::is_convertible_without_loss_to_integer(x[i])) {
      cpp11::stop("All values of '%s' must be integer-like, but '%s[%d]' was not",
                  name, name, i + 1);
    }
    ret[i] = static_cast<int>(std::round(x[i]));
  }
  return ret;
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
                                           bool require_later,
                                           const char * name) {
  auto time = to_vector_real<real_type>(r_time, name);
  auto prev = time_start;
  const auto eps = 1e-8;
  for (size_t i = 0; i < time.size(); ++i) {
    const auto t = time[i];
    if (!is_integer_like(t, eps)) {
      cpp11::stop("Expected 'time[%d]' to be integer-like", i + 1);
    }
    if (t < prev || (require_later && t == prev)) {
      cpp11::stop("Expected 'time[%d]' (%d) to be larger than the previous value (%d)",
                  i + 1, static_cast<int>(prev), static_cast<int>(t));
    }
    prev = t;
  }
  return time;
}

inline std::vector<size_t> check_index(cpp11::sexp r_index, size_t max,
                                       const char * name) {
  std::vector<size_t> ret;
  if (r_index == R_NilValue) {
    return ret;
  }
  if (TYPEOF(r_index) == REALSXP) {
    cpp11::doubles r_index_real = cpp11::as_cpp<cpp11::doubles>(r_index);
    return check_index(as_integers(r_index_real, name), max, name);
  }
  if (TYPEOF(r_index) != INTSXP) {
    cpp11::stop("Expected an integer vector for '%s'", name);
  }
  const int len = LENGTH(r_index);
  if (len == 0) {
    cpp11::stop("'%s' must have nonzero length", name);
  }
  ret.reserve(len);
  const auto data = INTEGER(r_index);
  for (int i = 0; i < len; ++i) {
    if (data[i] < 1 || data[i] > static_cast<int>(max)) {
      cpp11::stop("All values of '%s' must be in [1, %d], but '%s[%d]' was %d",
                  name, max, name, i + 1, data[i]);
    }
    ret.push_back(data[i] - 1);
  }
  return ret;
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
      const auto size_i = T::size_state(shared[i]);
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
void update_pars(T& obj, cpp11::list r_pars, bool grouped) {
  using system_type = typename T::system_type;
  if (grouped) {
    const auto n_groups = obj.n_groups();
    if (r_pars.size() != static_cast<int>(n_groups)) {
      cpp11::stop("Expected 'pars' to have length %d to match 'n_groups'",
                  static_cast<int>(n_groups));
    }
    for (size_t i = 0; i < n_groups; ++i) {
      cpp11::list r_pars_i = r_pars[i];
      obj.update_shared(i, [&] (auto& shared) {
                             system_type::update_shared(r_pars_i, shared);
                            });
    }
  } else {
    obj.update_shared(0, [&] (auto& shared) {
                           system_type::update_shared(r_pars, shared);
                          });
  }
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
        auto r_data_ij = cpp11::as_cpp<cpp11::list>(r_data_i[j]);
        data.push_back(T::build_data(r_data_ij));
      }
    }
  } else {
    for (size_t i = 0; i < n_time; ++i) {
      auto r_data_i = cpp11::as_cpp<cpp11::list>(r_data[i]);
      data.push_back(T::build_data(r_data_i));
    }
  }

  return data;
}

template <typename T>
void set_state(T& obj, cpp11::sexp r_state, bool grouped) {
  // Suppose that we have a n_state x n_particles x n_groups grouped
  // system, we then require that we have a state array with rank 3;
  // for an ungrouped system this will be rank 2 array.
  auto dim = cpp11::as_cpp<cpp11::integers>(r_state.attr("dim"));
  const auto rank = dim.size();
  const auto rank_expected = grouped ? 3 : 2;
  if (rank != rank_expected) {
    cpp11::stop("Expected 'state' to be a %dd array", rank_expected);
  }
  const int n_state = obj.n_state();
  const int n_particles =
    grouped ? obj.n_particles() : obj.n_particles() * obj.n_groups();
  const int n_groups = grouped ? obj.n_groups() : 1;
  if (dim[0] != n_state) {
    cpp11::stop("Expected the first dimension of 'state' to have size %d",
                n_state);
  }
  const auto recycle_particle = n_particles > 1 && dim[1] == 1;
  if (dim[1] != n_particles && dim[1] != 1) {
    cpp11::stop("Expected the second dimension of 'state' to have size %d or 1",
                n_particles);
  }

  const auto recycle_group = !grouped || (n_groups > 1 && dim[2] == 1);
  if (grouped && dim[2] != n_groups && dim[2] != 1) {
    cpp11::stop("Expected the third dimension of 'state' to have size %d or 1",
                n_groups);
  }
  obj.set_state(REAL(r_state), recycle_particle, recycle_group);
}

template <typename T>
SEXP rng_state_as_raw(const std::vector<T>& state) {
  const auto len = sizeof(T) * state.size();
  cpp11::writable::raws ret(len);
  std::memcpy(RAW(ret), state.data(), len);
  return ret;
}

// TODO: this will change later if we make state setting a bit more
// efficient.  Currently we accept a raw vector and copy that into an
// integer vetor, then copy that into chunks into the random number
// state. However, working with the direct raw vector is annoying, and
// we'll need to do some sort of reinterpret_cast on the pointer that
// RAW() returns I think?
template <typename rng_state_type>
inline auto check_rng_state(cpp11::sexp r_rng_state,
                            size_t n_streams,
                            const char * name) {
  using int_type = typename rng_state_type::int_type;
  const auto len = rng_state_type::size() * n_streams;
  const auto len_bytes = len * sizeof(int_type);

  if (!TYPEOF(r_rng_state)) {
    cpp11::stop("Expected a raw vector for '%s'", name);
  }
  cpp11::raws rng_state = cpp11::as_cpp<cpp11::raws>(r_rng_state);
  if (static_cast<size_t>(rng_state.size()) != len_bytes) {
    cpp11::stop("Incorrect length for '%s'; expected %d but given %d",
                name, static_cast<int>(len_bytes), rng_state.size());
  }
  std::vector<int_type> state(len);
  std::memcpy(state.data(), RAW(rng_state), len_bytes);
  return state;
}

inline double read_real(cpp11::list args, const char * name) {
  cpp11::sexp value = args[name];
  if (value == R_NilValue) {
    cpp11::stop("A value is expected for '%s'", name);
  }
  return to_double(value, name);
}

inline double read_real(cpp11::list args, const char * name,
                        double default_value) {
  cpp11::sexp value = args[name];
  return value == R_NilValue ? default_value : to_double(value, name);
}

inline int read_int(cpp11::list args, const char * name) {
  cpp11::sexp value = args[name];
  if (value == R_NilValue) {
    cpp11::stop("A value is expected for '%s'", name);
  }
  return to_int(value, name);
}

inline int read_int(cpp11::list args, const char * name,
                    int default_value) {
  cpp11::sexp value = args[name];
  return value == R_NilValue ? default_value : to_int(value, name);
}

inline size_t read_size(cpp11::list args, const char * name) {
  cpp11::sexp value = args[name];
  if (value == R_NilValue) {
    cpp11::stop("A value is expected for '%s'", name);
  }
  return to_size(value, name);
}

inline size_t read_size(cpp11::list args, const char * name,
                        size_t default_value) {
  cpp11::sexp value = args[name];
  return value == R_NilValue ? default_value : to_size(value, name);
}

inline bool read_bool(cpp11::list args, const char * name) {
  cpp11::sexp value = args[name];
  if (value == R_NilValue) {
    cpp11::stop("A value is expected for '%s'", name);
  }
  return to_bool(value, name);
}

inline bool read_bool(cpp11::list args, const char * name,
                      bool default_value) {
  cpp11::sexp value = args[name];
  return value == R_NilValue ? default_value : to_bool(value, name);
}

}
}
