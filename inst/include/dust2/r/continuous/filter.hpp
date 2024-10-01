#pragma once

#include <monty/r/random.hpp>
#include <dust2/r/helpers.hpp>
#include <dust2/r/filter.hpp>
#include <dust2/r/continuous/control.hpp>
#include <dust2/filter.hpp>
#include <dust2/continuous/system.hpp>

namespace dust2 {
namespace r {

template <typename T>
cpp11::sexp dust2_continuous_filter_alloc(cpp11::list r_pars,
                                          cpp11::sexp r_time_start,
                                          cpp11::sexp r_time,
                                          cpp11::list r_time_control,
                                          cpp11::list r_data,
                                          cpp11::sexp r_n_particles,
                                          cpp11::sexp r_n_groups,
                                          cpp11::sexp r_n_threads,
                                          cpp11::sexp r_index_state,
                                          cpp11::sexp r_seed) {
  using rng_state_type = typename T::rng_state_type;
  using rng_seed_type = std::vector<typename rng_state_type::int_type>;
  using real_type = typename  T::real_type;

  const auto n_particles = to_size(r_n_particles, "n_particles");
  const auto n_groups = to_size(r_n_groups, "n_groups");
  const auto n_threads = to_size(r_n_threads, "n_threads");
  const auto time_start = check_time(r_time_start, "time_start");
  const auto time = check_time_sequence(time_start, r_time, true, "time");
  const auto dt = check_dt(r_time_control, T::mixed_time::value);
  const auto ode_control = validate_ode_control<real_type>(r_time_control);
  const auto shared = build_shared<T>(r_pars, n_groups);
  const auto internal = build_internal<T>(shared, n_threads);
  const auto data = check_data<T>(r_data, shared, time.size(), "data");

  auto seed = monty::random::r::as_rng_seed<rng_state_type>(r_seed);
  const auto deterministic = false;

  // Create all the required rng states across the filter and the
  // system, in a reasonable way.  We need to make this slightly easier
  // to do from monty really.  Expand the state to give all the
  // state required by the filter (n_groups streams worth) and the
  // system (n_groups * n_particles worth, though the last bit of
  // expansion could be done by the system itself instead?)
  //
  // There are two ways of sorting out the state here:
  //
  // 1. we could take the first n_groups states for the filter and the
  // remaining for the systems.  This has the nice property that we can
  // expand the system state later if we support growing systems
  // (mrc-5355).  However, it has the undesirable consequence that a
  // filter with multiple groups will stream differently to a filter
  // containing a the first group only.
  //
  // 2. we take each block of (1+n_particles) states for each group,
  // giving the first to the filter and the rest to the system.  This
  // means that we can change the number of groups without affecting
  // the results, though we can't change the number of particles as
  // easily.
  const auto n_streams = n_groups * (n_particles + 1);
  const auto rng_state = monty::random::prng<rng_state_type>(n_streams, seed, deterministic).export_state();
  const auto rng_len = rng_state_type::size();
  rng_seed_type seed_filter;
  rng_seed_type seed_system;
  for (size_t i = 0; i < n_groups; ++i) {
    const auto it = rng_state.begin() + i * rng_len * (n_particles + 1);
    seed_filter.insert(seed_filter.end(),
                       it, it + rng_len);
    seed_system.insert(seed_system.end(),
                       it + rng_len, it + rng_len * (n_particles + 1));
  }

  const auto system = dust2::dust_continuous<T>(shared, internal, time_start, dt, ode_control, n_particles,
                                                seed_system, deterministic, n_threads);

  const auto index_state = check_index(r_index_state, system.n_state(),
                                       "index_state");

  auto obj = new filter<dust_continuous<T>>(system, time_start, time, data, index_state, seed_filter);
  cpp11::external_pointer<filter<dust_continuous<T>>> ptr(obj, true, false);

  cpp11::sexp r_n_state = cpp11::as_sexp(obj->sys.n_state());
  cpp11::sexp r_packing_state = packing_to_r(obj->sys.packing_state());

  using namespace cpp11::literals;
  return cpp11::writable::list{"ptr"_nm = ptr, "n_state"_nm = r_n_state, "packing_state"_nm = r_packing_state};
}

}
}
