#pragma once

#include <monty/r/random.hpp>
#include <dust2/r/helpers.hpp>
#include <dust2/r/unfilter.hpp>
#include <dust2/r/continuous/control.hpp>
#include <dust2/unfilter.hpp>
#include <dust2/continuous/system.hpp>

namespace dust2 {
namespace r {

// TODO: this name must be changed!
template <typename T>
cpp11::sexp dust2_continuous_unfilter_alloc(cpp11::list r_pars,
                                            cpp11::sexp r_time_start,
                                            cpp11::sexp r_time,
                                            cpp11::list r_time_control,
                                            cpp11::list r_data,
                                            cpp11::sexp r_n_particles,
                                            cpp11::sexp r_n_groups,
                                            cpp11::sexp r_n_threads) {
  using rng_state_type = typename T::rng_state_type;
  using real_type = typename  T::real_type;

  const auto n_particles = to_size(r_n_particles, "n_particles");
  const auto n_groups = to_size(r_n_groups, "n_groups");
  const auto n_threads = to_size(r_n_threads, "n_threads");
  const auto time_start = check_time(r_time_start, "time_start");
  const auto time = check_time_sequence(time_start, r_time, true, "time");
  const auto dt = check_dt(r_time_control, dust2::properties<T>::is_mixed_time::value, false);
  const auto ode_control = validate_ode_control<real_type>(r_time_control);
  const auto shared = build_shared<T>(r_pars, n_groups);
  const auto internal = build_internal<T>(shared, n_threads);
  const auto data = check_data<T>(r_data, shared, time.size(), "data");

  // It's possible that we don't want to always really be
  // deterministic here?  Though nooone can think of a case where
  // that's actually the behaviour wanted.  For now let's go fully
  // deterministic.
  auto seed = monty::random::seed_data<rng_state_type>(42);
  const auto deterministic = true;

  // Then allocate the system; this pulls together almost all the data
  // we need.  At this point we could have constructed the system out
  // of one that exists already on the R side, but I think that's
  // going to feel weirder overall.
  const auto system = dust2::dust_continuous<T>(shared, internal, time_start, dt, ode_control, n_particles,
                                                seed, deterministic, n_threads);

  auto obj = new unfilter<dust_continuous<T>>(system, time_start, time, data);
  cpp11::external_pointer<unfilter<dust_continuous<T>>> ptr(obj, true, false);

  cpp11::sexp r_n_state = cpp11::as_sexp(obj->sys.n_state());
  cpp11::sexp r_packing_state = packing_to_r(obj->sys.packing_state());
  cpp11::sexp r_packing_gradient = packing_to_r(obj->sys.packing_gradient());

  using namespace cpp11::literals;
  return cpp11::writable::list{"ptr"_nm = ptr, "n_state"_nm = r_n_state, "packing_state"_nm = r_packing_state, "packing_gradient"_nm = r_packing_gradient};
}

}
}
