#pragma once

#include <dust2/filter.hpp>
#include <mcstate/r/random.hpp>
#include <dust2/r/helpers.hpp>

namespace dust2 {
namespace r {

template <typename T>
cpp11::sexp dust2_filter_update_pars(cpp11::sexp ptr,
                                     cpp11::list r_pars) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  update_pars(obj->sys, r_pars);
  return R_NilValue;
}

template <typename T>
cpp11::sexp dust2_filter_run(cpp11::sexp ptr, cpp11::sexp r_initial,
                             bool save_history) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  if (r_initial != R_NilValue) {
    set_state(obj->sys, r_initial);
  }
  obj->run(r_initial == R_NilValue, save_history);

  cpp11::writable::doubles ret(obj->sys.n_groups());
  obj->last_log_likelihood(REAL(ret));
  return ret;
}

// Can collapse with above
template <typename T>
cpp11::sexp dust2_filter_last_history(cpp11::sexp ptr,
				      bool preserve_group_dimension) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  if (!obj->last_history_is_current()) {
    cpp11::stop("History is not current");
  }
  // We might relax this later, but will require some tools to work
  // with the output, really.
  constexpr bool reorder = true;

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
  if (preserve_group_dimension) {
    set_array_dims(ret, {n_state, n_particles, n_groups, n_times});
  } else {
    set_array_dims(ret, {n_state, n_particles * n_groups, n_times});
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_filter_rng_state(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  using rng_state_type = typename T::rng_state_type;

  // Undo the construction as above so that the rng state comes out in
  // the same format it goes in, as a single raw vector.
  const auto& state_filter = obj->rng_state();
  const auto& state_system = obj->sys.rng_state();
  const auto n_particles = obj->sys.n_particles();
  const auto n_groups = obj->sys.n_groups();
  const auto n_state = rng_state_type::size();
  const auto n_bytes = sizeof(typename rng_state_type::int_type);
  const auto n_bytes_state = n_bytes * n_state;
  cpp11::writable::raws ret(n_bytes * (state_filter.size() + state_system.size()));
  for (size_t i = 0; i < n_groups; ++i) {
    std::memcpy(RAW(ret) + i * n_bytes_state * (n_particles + 1),
                state_filter.data() + i * n_state,
                n_bytes_state);
    std::memcpy(RAW(ret) + i * n_bytes_state * (n_particles + 1) + n_bytes_state,
                state_system.data() + i * n_state * n_particles,
                n_bytes_state * n_particles);
  }

  return ret;
}

template <typename T>
cpp11::sexp dust2_filter_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  using rng_state_type = typename T::rng_state_type;
  using rng_int_type = typename rng_state_type::int_type;

  const auto n_particles = obj->sys.n_particles();
  const auto n_groups = obj->sys.n_groups();
  const auto n_streams = n_groups * (1 + n_particles);
  const auto n_state = rng_state_type::size();
  const auto rng_state =
    check_rng_state<rng_state_type>(r_rng_state, n_streams, "rng_state");

  std::vector<rng_int_type> state_filter(n_groups * n_state);
  std::vector<rng_int_type> state_system(n_groups * n_particles * n_state);
  for (size_t i = 0; i < n_groups; ++i) {
    const auto src = rng_state.begin() + i * (1 + n_particles) * n_state;
    std::copy_n(src,
                n_state,
                state_filter.begin() + i * n_state);
    std::copy_n(src + n_state,
                n_state * n_particles,
                state_system.begin() + i * n_state * n_particles);
  }

  obj->set_rng_state(state_filter);
  obj->sys.set_rng_state(state_system);

  return R_NilValue;
}

}
}
