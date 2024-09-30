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

  auto ode_control = validate_ode_control<real_type>(r_time_control);

  auto obj = new dust_continuous<T>(shared, internal, time, ode_control,
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
                                  bool include_coefficients) {
  using namespace cpp11::literals;
  auto ret = cpp11::writable::list{
    "dydt"_nm = cpp11::as_sexp(internals.dydt),
    "step_times"_nm = cpp11::as_sexp(internals.step_times),
    "step_size"_nm = cpp11::as_sexp(internals.step_size),
    "error"_nm = cpp11::as_sexp(internals.error),
    "n_steps"_nm = cpp11::as_sexp(internals.n_steps),
    "n_steps_accepted"_nm = cpp11::as_sexp(internals.n_steps_accepted),
    "n_steps_rejected"_nm = cpp11::as_sexp(internals.n_steps_rejected),
    "coefficients"_nm = R_NilValue};
  if (include_coefficients) {
    auto r_coef = cpp11::writable::doubles_matrix<>(internals.c1.size(), 5);
    auto coef = REAL(r_coef);
    coef = std::copy(internals.c1.begin(), internals.c1.end(), coef);
    coef = std::copy(internals.c2.begin(), internals.c2.end(), coef);
    coef = std::copy(internals.c3.begin(), internals.c3.end(), coef);
    coef = std::copy(internals.c4.begin(), internals.c4.end(), coef);
    coef = std::copy(internals.c5.begin(), internals.c5.end(), coef);
    ret["coefficients"] = r_coef;
  }
  return ret;
}

template <typename T>
SEXP dust2_system_internals(cpp11::sexp ptr, bool include_coefficients) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const auto& internals = obj->ode_internals();
  cpp11::writable::list ret(internals.size());
  for (size_t i = 0; i < internals.size(); ++i) {
    ret[i] = ode_internals_to_sexp(internals[i], include_coefficients);
  }
  return cpp11::as_sexp(ret);
}

}
}
