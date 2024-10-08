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
    cpp11::stop("Expected '%s' to have length %d, but it had length %d",
                name, len, LENGTH(x));
  }
}

inline void check_rank(cpp11::sexp value, size_t rank, const char * name) {
  auto r_dim = value.attr("dim");
  size_t rank_given = r_dim == R_NilValue ? 1 : LENGTH(r_dim);

  if (rank_given == rank) {
    return;
  }

  if (rank == 1) {
    if (rank_given == 2) {
      cpp11::stop("Expected '%s' to be a vector, but was given a matrix", name);
    } else {
      cpp11::stop("Expected '%s' to be a vector, but was given a %d-dimensional array",
                  name, static_cast<int>(rank_given));
    }
  } else if (rank == 2) {
    if (rank_given == 1) {
      cpp11::stop("Expected '%s' to be a matrix, but was given a vector", name);
    } else {
      cpp11::stop("Expected '%s' to be a matrix, but was given a %d-dimensional array",
                  name, static_cast<int>(rank_given));
    }
  } else {
    if (rank_given == 1) {
      cpp11::stop("Expected '%s' to be a %d-dimensional array, but was given a vector",
                  name, static_cast<int>(rank));
    } else if (rank_given == 2) {
      cpp11::stop("Expected '%s' to be a %d-dimensional array, but was given a matrix",
                  name, static_cast<int>(rank));
    } else {
      cpp11::stop("Expected '%s' to be a %d-dimensional array, but was given a %d-dimensional array",
                  name, static_cast<int>(rank), static_cast<int>(rank_given));
    }
  }
}

template <size_t rank>
void check_dimensions(cpp11::sexp value,
                      const dust2::array::dimensions<rank>& dim,
                      const char * name) {
  check_rank(value, rank, name);
  if (rank == 1) {
    check_length(value, dim.size, name);
  } else {
    auto r_dim = cpp11::as_cpp<cpp11::integers>(value.attr("dim"));
    for (size_t i = 0; i < rank; ++i) {
      if (static_cast<size_t>(r_dim[i]) != dim.dim[i]) {
        if (rank > 2) {
          cpp11::stop("Expected dimension %d of '%s' to have length %d, but it had length %d",
                      static_cast<int>(i + 1), name, static_cast<int>(dim.dim[i]), r_dim[i]);
        } else if (i == 0) {
          cpp11::stop("Expected '%s' to have %d rows, but was given %d",
                      name, static_cast<int>(dim.dim[i]), r_dim[i]);
        } else {
          cpp11::stop("Expected '%s' to have %d columns, but was given %d",
                      name, static_cast<int>(dim.dim[i]), r_dim[i]);
        }
      }
    }
  }
}


template <size_t rank>
dust2::array::dimensions<rank> read_dimensions(cpp11::list args,
                                               const char * name) {
  cpp11::sexp value = args[name];
  check_rank(value, rank, name);
  auto r_dim = cpp11::as_cpp<cpp11::integers>(value.attr("dim"));
  return dust2::array::dimensions<rank>(r_dim.begin());
}

template <>
inline dust2::array::dimensions<1> read_dimensions(cpp11::list args,
                                                   const char * name) {
  SEXP value = args[name];
  check_rank(value, 1, name);
  const size_t len = LENGTH(value);
  return dust2::array::dimensions<1>{len};
}

inline double to_double(cpp11::sexp x, bool allow_na, const char * name) {
  check_scalar(x, name);
  if (TYPEOF(x) == REALSXP) {
    double ret = REAL(x)[0];
    if (!allow_na && ISNA(ret)) {
      cpp11::stop("'%s' must not be a missing value (NA)", name);
    }
    return ret;
  }
  if (TYPEOF(x) == INTSXP) {
    int ret = INTEGER(x)[0];
    if (ISNA(ret)) {
      if (allow_na) {
        return NA_REAL;
      } else {
        cpp11::stop("'%s' must not be a missing value (NA)", name);
      }
    }
    return static_cast<double>(ret);
  }
  if (allow_na && TYPEOF(x) == LGLSXP && ISNA(INTEGER(x)[0])) {
    return NA_REAL;
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
  const auto allow_na = false;
  return to_double(r_time, allow_na, name);
}


inline double check_dt(cpp11::list r_time_control, bool enabled = true,
                       bool required = true) {
  cpp11::sexp r_dt = r_time_control["dt"];
  const bool is_null = r_dt == R_NilValue;

  if (is_null && !required) {
    return 0;
  }
  if (!is_null && !enabled) {
    cpp11::stop("Unexpected value for 'dt' in time control");
  }

  const auto allow_na = false;
  const auto dt = to_double(r_dt, allow_na, "dt");
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
  if (time[0] < time_start || (require_later && time[0] == time_start)) {
    cpp11::stop("Expected 'time[1]' (%d) to be larger than the previous value (%d)",
                static_cast<int>(time[0]), static_cast<int>(time_start));
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
  if (dims.size() < 2) {
    return;
  }
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
  if (r_pars.size() != static_cast<int>(n_groups)) {
    cpp11::stop("Expected 'pars' to have length %d to match 'n_groups'",
                static_cast<int>(n_groups));
  }

  std::vector<typename T::shared_state> shared;
  shared.reserve(n_groups);

  // TODO: we could check gradient here too but it seems a bit over
  // the top.
  dust2::packing packing{};
  for (size_t i = 0; i < n_groups; ++i) {
    shared.push_back(T::build_shared(r_pars[i]));
    const auto packing_i = T::packing_state(shared[i]);
    if (i == 0) {
      packing = packing_i;
    } else if (packing_i != packing) {
      const auto size_i = packing_i.size();
      const auto size = packing.size();
      // Later, we might do better with these errors, but practically
      // this means throwing an R error from C++ so that we can get
      // this data there to do the comparison with and build strings
      // from.
      if (size_i != size) {
        cpp11::stop("State length for group %d was different to previous groups; total length was expected to be %d but it was %d",
                    i + 1, size, size_i);
      } else {
        cpp11::stop("State length for group %d was different to previous groups; total length %d as expected, but the internal structure was different",
                    i + 1, size);
      }
    }
  }

  return shared;
}

template <typename T>
std::vector<typename T::internal_state> build_internal(const std::vector<typename T::shared_state>& shared, size_t n_threads) {
  std::vector<typename T::internal_state> internal;
  const size_t n_groups = shared.size();
  internal.reserve(n_groups * n_threads);
  // We can parallelise this but it's probably not really wanted.
  for (size_t i = 0; i < n_threads; ++i) {
    for (auto& s : shared) {
      internal.push_back(T::build_internal(s));
    }
  }
  return internal;
}

template <typename T>
void update_pars(T& obj, cpp11::list r_pars, const std::vector<size_t>& index_group) {
  using system_type = typename T::system_type;

  const auto n_groups = index_group.size();
  if (r_pars.size() != static_cast<int>(n_groups)) {
    cpp11::stop("Expected 'pars' to have length %d to match '%s'",
                static_cast<int>(n_groups),
                index_group.size() == obj.n_groups() ? "n_groups" : "index_group");
  }
  for (size_t i = 0; i < n_groups; ++i) {
    cpp11::list r_pars_i = r_pars[i];
    obj.update_shared(index_group[i], [&] (auto& shared) {
      system_type::update_shared(r_pars_i, shared);
    });
  }
}

template <typename T>
std::vector<typename T::data_type> check_data(cpp11::list r_data,
                                              const std::vector<typename T::shared_state>& shared,
                                              size_t n_time,
                                              const char * name) {
  const size_t n_groups = shared.size();
  std::vector<typename T::data_type> data;
  data.reserve(n_time * n_groups);

  // Errors here are no longer likely to be thrown, as the R-side
  // check_data() should do the work for us.  The exception is that
  // T::build_data() might fail and we should probably report back the
  // time and group index that is failing here.
  check_length(r_data, n_time, name);
  for (size_t i = 0; i < n_time; ++i) {
    auto r_data_i = cpp11::as_cpp<cpp11::list>(r_data[i]);
    check_length(r_data_i, n_groups, "data[i]");
    for (size_t j = 0; j < n_groups; ++j) {
      auto r_data_ij = cpp11::as_cpp<cpp11::list>(r_data_i[j]);
      data.push_back(T::build_data(r_data_ij, shared[j]));
    }
  }

  return data;
}

template <typename T>
void set_state(T& obj, cpp11::sexp r_state, bool preserve_group_dimension,
               const std::vector<size_t>& index_group) {
  // Suppose that we have a n_state x n_particles x n_groups grouped
  // system, we then require that we have a state array with rank 3;
  // for an ungrouped system this will be rank 2 array.
  //
  // TODO: these checks would be nicer in R, just do it there and here
  // we can just accept what we are given? (mrc-5565)
  auto dim = cpp11::as_cpp<cpp11::integers>(r_state.attr("dim"));
  const auto rank = dim.size();
  const auto rank_expected = preserve_group_dimension ? 3 : 2;
  if (rank != rank_expected) {
    cpp11::stop("Expected 'state' to be a %dd array", rank_expected);
  }
  const int n_state = obj.n_state();
  const int n_particles =
    preserve_group_dimension ? obj.n_particles() :
    obj.n_particles() * obj.n_groups();
  const int n_groups = index_group.size();
  if (dim[0] != n_state) {
    cpp11::stop("Expected the first dimension of 'state' to have size %d",
                n_state);
  }
  const auto recycle_particle = n_particles > 1 && dim[1] == 1;
  if (dim[1] != n_particles && dim[1] != 1) {
    cpp11::stop("Expected the second dimension of 'state' to have size %d or 1",
                n_particles);
  }

  const auto recycle_group =
    !preserve_group_dimension || (n_groups > 1 && dim[2] == 1);
  if (preserve_group_dimension && dim[2] != n_groups && dim[2] != 1) {
    cpp11::stop("Expected the third dimension of 'state' to have size %d or 1",
                n_groups);
  }
  obj.set_state(REAL(r_state), recycle_particle, recycle_group,
                index_group);
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
  const bool allow_na = false;
  cpp11::sexp value = args[name];
  if (value == R_NilValue) {
    cpp11::stop("A value is expected for '%s'", name);
  }
  return to_double(value, allow_na, name);
}

inline double read_real(cpp11::list args, const char * name,
                        double default_value) {
  const bool allow_na = ISNA(default_value);
  cpp11::sexp value = args[name];
  return value == R_NilValue ? default_value : to_double(value, allow_na, name);
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

template <typename real_type>
void read_real_vector(cpp11::list args, size_t len, real_type * dest,
                      const char *name, bool required) {
  cpp11::sexp value = args[name];
  if (value == R_NilValue) {
    if (required) {
      cpp11::stop("A value is expected for '%s'", name);
    }
  } else {
    check_length(value, len, name);
    if (TYPEOF(value) == REALSXP) {
      auto value_dbl = cpp11::as_cpp<cpp11::doubles>(value);
      std::copy(value_dbl.begin(), value_dbl.end(), dest);
    } else if (TYPEOF(value) == INTSXP) {
      auto value_int = cpp11::as_cpp<cpp11::integers>(value);
      std::copy(value_int.begin(), value_int.end(), dest);
    } else {
      cpp11::stop("'%s' must be numeric", name);
    }
  }
}

template <typename real_type, size_t rank>
void read_real_array(cpp11::list args, const dust2::array::dimensions<rank>& dim,
                     real_type * dest, const char *name, bool required) {
  cpp11::sexp value = args[name];
  if (value == R_NilValue) {
    if (required) {
      cpp11::stop("A value is expected for '%s'", name);
    }
  } else {
    check_dimensions<rank>(value, dim, name);
    if (TYPEOF(value) == REALSXP) {
      auto value_dbl = cpp11::as_cpp<cpp11::doubles>(value);
      std::copy(value_dbl.begin(), value_dbl.end(), dest);
    } else if (TYPEOF(value) == INTSXP) {
      auto value_int = cpp11::as_cpp<cpp11::integers>(value);
      std::copy(value_int.begin(), value_int.end(), dest);
    } else {
      cpp11::stop("'%s' must be numeric", name);
    }
  }
}

inline cpp11::sexp packing_to_r(const dust2::packing& packing_info) {
  if (packing_info.size() == 0) {
    return R_NilValue;
  }
  using namespace cpp11::literals;

  const auto& data = packing_info.data();
  const auto len = data.size();
  cpp11::writable::list ret(len);
  cpp11::writable::strings nms(len);
  for (size_t i = 0; i < len; ++i) {
    nms[i] = data[i].first;
    ret[i] = cpp11::writable::integers(data[i].second.begin(),
                                       data[i].second.end());
  }
  ret.attr("names") = nms;
  return ret;
}

}
}
