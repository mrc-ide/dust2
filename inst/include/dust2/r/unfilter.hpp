#pragma once

#include <dust2/unfilter.hpp>
#include <mcstate/r/random.hpp>
#include <dust2/r/helpers.hpp>

namespace dust2 {
namespace r {

template <typename T>
cpp11::sexp dust2_unfilter_update_pars(cpp11::sexp ptr,
                                       cpp11::list r_pars,
                                       bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter<T>>>(ptr).get();
  update_pars(obj->sys, cpp11::as_cpp<cpp11::list>(r_pars), grouped);
  return R_NilValue;
}

template <typename T>
cpp11::sexp dust2_unfilter_run(cpp11::sexp ptr, cpp11::sexp r_initial,
                               bool save_history, bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter<T>>>(ptr).get();
  if (r_initial != R_NilValue) {
    set_state(obj->sys, r_initial, grouped);
  }
  obj->run(r_initial == R_NilValue, save_history);

  const auto n_groups = obj->sys.n_groups();
  const auto n_particles = obj->sys.n_particles();
  cpp11::writable::doubles ret(n_groups * n_particles);
  obj->last_log_likelihood(REAL(ret));
  if (grouped && n_particles > 1) {
    set_array_dims(ret, {n_particles, n_groups});
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_unfilter_last_history(cpp11::sexp ptr, bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter<T>>>(ptr).get();
  if (!obj->last_history_is_current()) {
    cpp11::stop("History is not current");
  }

  constexpr bool reorder = false; // never needed

  const auto& history = obj->last_history();
  const auto& dims = history.dims();
  // Could use destructured bind here in recent C++?
  const auto n_state = dims[0];
  const auto n_particles = dims[1];
  const auto n_groups = dims[2];
  const auto n_times = dims[3];
  const auto len = n_state * n_particles * n_groups * n_times;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  history.export_state(REAL(ret), reorder);
  if (grouped) {
    set_array_dims(ret, {n_state, n_particles, n_groups, n_times});
  } else {
    set_array_dims(ret, {n_state, n_particles * n_groups, n_times});
  }
  return ret;
}

}
}
