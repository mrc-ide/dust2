#pragma once

namespace dust2 {
namespace r {

inline void check_scalar(cpp11::sexp x, const char * name) {
  if (LENGTH(x) != 1) {
    cpp11::stop("'%s' must be scalar", name);
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
