#pragma once

namespace dust2 {
namespace r {

inline void check_scalar(cpp11::sexp x, const char * name) {
  if (LENGTH(x) != 1) {
    cpp11::stop("'%s' must be a scalar", name);
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
  return std::abs(x - round(x)) <= eps;
}

inline double check_time(cpp11::sexp r_time) {
  const auto time = to_double(r_time, "time");
  const auto eps = 1e-8;
  // We can relax this later and carefully align time onto a grid
  if (!is_integer_like(time, eps)) {
    throw std::runtime_error("Expected 'time' to be integer-like");
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
    throw std::runtime_error("Expected 'dt' to be the inverse of an integer");
  }
  return dt;
}

// template <typename real_type>
// inline cpp11::sexp to_matrix(std::vector<real_type> x, size_t nr, size_t nc) {
//   cpp11::writable::integers dim{static_cast<int>(nr), static_cast<int>(nc)};
//   cpp11::writable::doubles ret(x.size());
//   std::copy(x.begin(), x.end(), REAL(ret));
//   ret.attr("dim") = dim;
//   return ret;
// }

}
}
