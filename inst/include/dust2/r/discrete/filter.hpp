#pragma once

#include <mcstate/r/random.hpp>
#include <dust2/r/helpers.hpp>
#include <dust2/r/filter.hpp>
#include <dust2/filter.hpp>
#include <dust2/discrete/system.hpp>

namespace dust2 {
namespace r {

template <typename T>
cpp11::sexp dust2_discrete_filter_alloc(cpp11::list r_pars,
                                        cpp11::sexp r_time_start,
                                        cpp11::sexp r_time,
                                        cpp11::sexp r_dt,
                                        cpp11::list r_data,
                                        cpp11::sexp r_n_particles,
                                        cpp11::sexp r_n_groups,
                                        cpp11::sexp r_index,
                                        cpp11::sexp r_seed) {
  using rng_state_type = typename T::rng_state_type;
  using rng_seed_type = std::vector<typename rng_state_type::int_type>;

  const auto n_particles = to_size(r_n_particles, "n_particles");
  const auto n_groups = to_size(r_n_groups, "n_groups");
  const auto grouped = n_groups > 0;
  const auto time_start = check_time(r_time_start, "time_start");
  const auto time = check_time_sequence(time_start, r_time, true, "time");
  const auto dt = check_dt(r_dt);
  const auto shared = build_shared<T>(r_pars, n_groups);
  const auto internal = build_internal<T>(shared);
  const auto data = check_data<T>(r_data, time.size(), n_groups, "data");

  // It's possible that we don't want to always really be
  // deterministic here?  Though nooone can think of a case where
  // that's actually the behaviour wanted.  For now let's go fully
  // deterministic.
  auto seed = mcstate::random::r::as_rng_seed<rng_state_type>(r_seed);
  const auto deterministic = false;

  // Create all the required rng states across the filter and the
  // system, in a reasonable way.  We need to make this slightly easier
  // to do from mcstate really.  Expand the state to give all the
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
  const auto n_groups_effective = grouped ? n_groups : 1;
  const auto n_streams = n_groups_effective * (n_particles + 1);
  const auto rng_state = mcstate::random::prng<rng_state_type>(n_streams, seed, deterministic).export_state();
  const auto rng_len = rng_state_type::size();
  rng_seed_type seed_filter;
  rng_seed_type seed_system;
  for (size_t i = 0; i < n_groups_effective; ++i) {
    const auto it = rng_state.begin() + i * rng_len * (n_particles + 1);
    seed_filter.insert(seed_filter.end(),
                       it, it + rng_len);
    seed_system.insert(seed_system.end(),
                      it + rng_len, it + rng_len * (n_particles + 1));
  }

  const size_t n_threads = 1;
  const auto system = dust2::dust_discrete<T>(shared, internal, time_start, dt, n_particles,
                                              seed_system, deterministic, n_threads);

  const auto index = check_index(r_index, system.n_state(), "index");

  auto obj = new filter<dust_discrete<T>>(system, time_start, time, data, index, seed_filter);
  cpp11::external_pointer<filter<dust_discrete<T>>> ptr(obj, true, false);

  cpp11::sexp r_n_state = cpp11::as_sexp(obj->sys.n_state());
  cpp11::sexp r_group_names = R_NilValue;
  if (grouped) {
    r_group_names = r_pars.attr("names");
  }
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
