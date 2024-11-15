#include <cpp11.hpp>

[[cpp11::register]]
SEXP check_externalptr(SEXP ptr, SEXP allow_rebuild) {
  if (ptr == R_NilValue) {
    return Rf_ScalarLogical(false);
  } else if (R_ExternalPtrAddr(ptr)) {
    return Rf_ScalarLogical(true);
  } else {
    if (!INTEGER(allow_rebuild)[0]) {
      cpp11::stop("Pointer has been serialised, cannot continue safely");
    }
    return Rf_ScalarLogical(false);
  }
}
