#pragma once

#include <dust2/r/helpers.hpp>
#include <monty/r/random.hpp>
#include <dust2/r/system.hpp>
#include <dust2/r/continuous/control.hpp>

#include <dust2/continuous/system.hpp>

namespace dust2 {
namespace r {

template <typename T>
SEXP dust2_continuous_alloc(cpp11::list r_pars,
                            cpp11::sexp r_time,
                            cpp11::list r_time_control,
                            cpp11::sexp r_n_particles,
                            cpp11::sexp r_n_groups,
                            cpp11::sexp r_seed,
                            cpp11::sexp r_deterministic,
                            cpp11::sexp r_n_threads) {
  using rng_state_type = typename T::rng_state_type;
  using real_type = typename  T::real_type;

  const auto time = check_time(r_time, "time");

  const auto n_particles = to_size(r_n_particles, "n_particles");
  const auto n_groups = to_size(r_n_groups, "n_groups");
  const auto n_threads = to_size(r_n_threads, "n_threads");

  const auto shared = build_shared<T>(r_pars, n_groups);
  const auto internal = build_internal<T>(shared, n_threads);

  auto seed = monty::random::r::as_rng_seed<rng_state_type>(r_seed);
  auto deterministic = to_bool(r_deterministic, "deterministic");

  auto dt = check_dt(r_time_control, dust2::properties<T>::is_mixed_time::value, false);
  auto ode_control = validate_ode_control<real_type>(r_time_control);

  auto obj = new dust_continuous<T>(shared, internal, time, dt, ode_control,
                                    n_particles, seed, deterministic,
                                    n_threads);
  cpp11::external_pointer<dust_continuous<T>> ptr(obj, true, false);

  // Later, we'll export a bit more back from the system (in particular
  // systems need to provide information about how they organise
  // variables.
  cpp11::sexp r_n_state = cpp11::as_sexp(obj->n_state());
  cpp11::sexp r_packing_state = packing_to_r(obj->packing_state());
  cpp11::sexp r_packing_gradient = packing_to_r(obj->packing_gradient());

  using namespace cpp11::literals;
  return cpp11::writable::list{"ptr"_nm = ptr, "n_state"_nm = r_n_state, "packing_state"_nm = r_packing_state, "packing_gradient"_nm = r_packing_gradient};
}

template <typename real_type>
cpp11::sexp ode_internals_to_sexp(const ode::internals<real_type>& internals,
                                  bool include_coefficients,
                                  bool include_history,
                                  std::vector<std::string> event_names) {
  using namespace cpp11::literals;
  auto ret = cpp11::writable::list{
    "dydt"_nm = cpp11::as_sexp(internals.dydt),
    "step_times"_nm = cpp11::as_sexp(internals.step_times),
    "step_size"_nm = cpp11::as_sexp(internals.step_size),
    "error"_nm = cpp11::as_sexp(internals.error),
    "n_steps"_nm = cpp11::as_sexp(internals.n_steps),
    "n_steps_accepted"_nm = cpp11::as_sexp(internals.n_steps_accepted),
    "n_steps_rejected"_nm = cpp11::as_sexp(internals.n_steps_rejected),
    "coefficients"_nm = R_NilValue,
    "history"_nm = R_NilValue,
    "events"_nm = R_NilValue};

  const auto n_variables = internals.last.c2.size();

  if (include_coefficients) {
    auto r_coef = cpp11::writable::doubles_matrix<>(internals.last.c2.size(), 5);
    auto coef = REAL(r_coef);
    // TODO: we might want to use copy_n everywhere, just a little
    // more churn.  Perhaps internals should also advertise
    // n_variables more obviously?
    coef = std::copy_n(internals.last.c1.begin(), n_variables, coef);
    coef = std::copy(internals.last.c2.begin(), internals.last.c2.end(), coef);
    coef = std::copy(internals.last.c3.begin(), internals.last.c3.end(), coef);
    coef = std::copy(internals.last.c4.begin(), internals.last.c4.end(), coef);
    coef = std::copy(internals.last.c5.begin(), internals.last.c5.end(), coef);
    ret["coefficients"] = r_coef;
  }
  if (include_history && !internals.history_values.empty()) {
    const auto n_history_entries = internals.history_values.size();
    const auto n_history_state = internals.history_values.n_state();
    auto r_history_coef =
      cpp11::writable::doubles(n_history_state * 5 * n_history_entries);
    auto r_history_t0 = cpp11::writable::doubles(n_history_entries);
    auto r_history_t1 = cpp11::writable::doubles(n_history_entries);
    auto r_history_h = cpp11::writable::doubles(n_history_entries);
    auto history_coef = REAL(r_history_coef);
    auto history_t0 = REAL(r_history_t0);
    auto history_t1 = REAL(r_history_t1);
    auto history_h = REAL(r_history_h);
    for (auto& h: internals.history_values.data()) {
      *history_t0++ = h.t0;
      *history_t1++ = h.t1;
      *history_h++ = h.h;
      history_coef = std::copy_n(h.c1.begin(), n_variables, history_coef);
      history_coef = std::copy(h.c2.begin(), h.c2.end(), history_coef);
      history_coef = std::copy(h.c3.begin(), h.c3.end(), history_coef);
      history_coef = std::copy(h.c4.begin(), h.c4.end(), history_coef);
      history_coef = std::copy(h.c5.begin(), h.c5.end(), history_coef);
    }
    set_array_dims(r_history_coef, {n_history_state, 5, n_history_entries});
    auto r_history = cpp11::writable::list{"t0"_nm = std::move(r_history_t0),
                                           "t1"_nm = std::move(r_history_t1),
                                           "h"_nm = std::move(r_history_h),
                                           "coefficients"_nm = std::move(r_history_coef)};
    ret["history"] = cpp11::as_sexp(r_history);
  }
  if (!event_names.empty()) {
    const auto n_events = internals.events.size();
    auto r_event_time = cpp11::writable::doubles(n_events);
    auto r_event_index = cpp11::writable::integers(n_events);
    auto r_event_sign = cpp11::writable::doubles(n_events);
    auto r_event_name = cpp11::writable::strings(n_events);
    for (size_t i = 0; i < n_events; ++i) {
      const auto idx = internals.events[i].index;
      r_event_time[i] = internals.events[i].time;
      r_event_index[i] = static_cast<int>(idx) + 1;
      r_event_sign[i] = internals.events[i].sign;
      r_event_name[i] = event_names[idx];
    }
    auto r_events = cpp11::writable::list{"name"_nm = std::move(r_event_name),
                                          "time"_nm = std::move(r_event_time),
                                          "index"_nm = std::move(r_event_index),
                                          "sign"_nm = std::move(r_event_sign)};
    ret["events"] = cpp11::as_sexp(r_events);
  }
  return ret;
}

template <typename T>
SEXP dust2_system_internals(cpp11::sexp ptr, bool include_coefficients, bool include_history) {
  auto *obj = dust2::r::safely_read_externalptr<T>(ptr, "system_internals");
  const auto& internals = obj->ode_internals();
  cpp11::writable::list ret(internals.size());
  const auto n_groups = obj->n_groups();
  const auto n_particles = obj->n_particles();
  for (size_t i = 0; i < n_groups; ++i) {
    const auto event_names = obj->event_names(i);
    for (size_t j = 0; j < n_particles; ++j) {
      const auto k = i * n_particles + j;
      ret[k] = ode_internals_to_sexp(internals[k],
                                     include_coefficients,
                                     include_history,
                                     event_names);
    }
  }
  return cpp11::as_sexp(ret);

}

}
}
