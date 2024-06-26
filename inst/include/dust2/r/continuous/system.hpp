#pragma once

#include <dust2/r/helpers.hpp>
#include <mcstate/r/random.hpp>
#include <dust2/r/system.hpp>
#include <dust2/r/continuous/control.hpp>

#include <dust2/continuous/system.hpp>

namespace dust2 {
namespace r {

template <typename T>
SEXP dust2_continuous_alloc(cpp11::list r_pars,
                            cpp11::sexp r_time,
                            cpp11::list r_ode_control,
                            cpp11::sexp r_n_particles,
                            cpp11::sexp r_n_groups,
                            cpp11::sexp r_seed,
                            cpp11::sexp r_deterministic) {
  using rng_state_type = typename T::rng_state_type;
  using real_type = typename  T::real_type;

  const auto time = check_time(r_time, "time");

  const auto n_particles = to_size(r_n_particles, "n_particles");
  const auto n_groups = to_size(r_n_groups, "n_groups");

  const auto shared = build_shared<T>(r_pars, n_groups);
  // Later, we need one of these per thread
  const auto internal = build_internal<T>(shared);

  auto seed = mcstate::random::r::as_rng_seed<rng_state_type>(r_seed);
  auto deterministic = to_bool(r_deterministic, "deterministic");
  auto ode_control = validate_ode_control<real_type>(r_ode_control);

  auto obj = new dust_continuous<T>(shared, internal, time, ode_control,
                                    n_particles, seed, deterministic);
  cpp11::external_pointer<dust_continuous<T>> ptr(obj, true, false);

  // Later, we'll export a bit more back from the system (in particular
  // systems need to provide information about how they organise
  // variables.
  const auto grouped = n_groups > 0;
  cpp11::sexp r_group_names = R_NilValue;
  if (grouped) {
    r_group_names = r_pars.attr("names");
  }
  cpp11::sexp r_n_state = cpp11::as_sexp(obj->n_state());
  cpp11::sexp r_grouped = cpp11::as_sexp(grouped);

  using namespace cpp11::literals;
  return cpp11::writable::list{"ptr"_nm = ptr,
      "n_state"_nm = r_n_state,
      "grouped"_nm = r_grouped,
      "group_names"_nm = r_group_names
      };
}

}
}
