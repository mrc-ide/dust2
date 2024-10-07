#pragma once

#include <cstring>
#include <dust2/common.hpp>
#include <cpp11.hpp>

namespace dust2 {
namespace r {

namespace debug {
inline cpp11::environment create_env() {
  cpp11::function new_env(cpp11::package("base")["new.env"]);
  return new_env();
}

template <typename T>
void save(T x, const char * name, cpp11::environment env) {
  env[name] = x;
}

template <typename T, size_t rank>
void save(T *x, array::dimensions<rank> dim, const char * name,
          cpp11::environment env) {
  const size_t size = dim.size();
  cpp11::writable::doubles r_x(size);
  std::copy_n(x, size, REAL(r_x));
  if (rank > 1) {
    cpp11::writable::integers r_dim(rank);
    std::copy_n(dim.dim.begin(), rank, INTEGER(r_dim));
    r_x.attr("dim") = r_dim;
  }
  env[name] = r_x;
}

inline void browser(cpp11::environment env, const char *phase, double time) {
  cpp11::function debug_env(cpp11::package("dust2")["debug_env"]);
  cpp11::sexp r_phase = cpp11::as_sexp(phase);
  cpp11::sexp r_time = cpp11::as_sexp(time);
  debug_env(env, r_phase, r_time);
}

}
}
}
