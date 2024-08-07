#pragma once

#include <dust2/r/helpers.hpp>
#include <dust2/r/system.hpp>
#include <mcstate/r/random.hpp>

#include <dust2/discrete/system.hpp>

namespace dust2 {
namespace r {

template <typename T>
SEXP dust2_discrete_alloc(cpp11::list r_pars,
                          cpp11::sexp r_time,
                          cpp11::sexp r_dt,
                          cpp11::sexp r_n_particles,
                          cpp11::sexp r_n_groups,
                          cpp11::sexp r_seed,
                          cpp11::sexp r_deterministic,
			  cpp11::sexp r_n_threads) {
  using rng_state_type = typename T::rng_state_type;

  // These duplicate checks that happen on the R side and can be
  // relaxed over time.  However, they're fast and pretty harmless -
  // most will just arrange for the conversion from SEXP type to the
  // expected underlying C type (with no cost) and then we cast out
  // into the C++ type that we need here (e.g., SEXP -> int -> size_t)
  //
  // The only one of these that might throw is the the build_shared
  // call, which might fail because of validation checks that the user
  // has getting parameters from the SEXP, or because it creates a
  // system with an unexpected size.  These will result in a fairly
  // ugly error compared with most.
  const auto time = check_time(r_time, "time");
  const auto dt = check_dt(r_dt);
  const auto n_particles = to_size(r_n_particles, "n_particles");
  const auto n_groups = to_size(r_n_groups, "n_groups");
  const auto n_threads = to_size(r_n_threads, "n_threads");
  const auto shared = build_shared<T>(r_pars, n_groups);
  const auto internal = build_internal<T>(shared, n_threads);
  auto seed = mcstate::random::r::as_rng_seed<rng_state_type>(r_seed);
  auto deterministic = to_bool(r_deterministic, "deterministic");

  auto obj = new dust_discrete<T>(shared, internal, time, dt, n_particles,
                                  seed, deterministic, n_threads);
  cpp11::external_pointer<dust_discrete<T>> ptr(obj, true, false);

  // Later, we'll export information about how systems structure
  // variables (mrc-5422, with support needed from mcstate2)
  cpp11::sexp r_n_state = cpp11::as_sexp(obj->n_state());

  using namespace cpp11::literals;
  return cpp11::writable::list{"ptr"_nm = ptr, "n_state"_nm = r_n_state};
}

}
}
