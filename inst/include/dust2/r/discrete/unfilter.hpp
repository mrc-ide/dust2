#pragma once

#include <mcstate/r/random.hpp>
#include <dust2/r/helpers.hpp>
#include <dust2/r/unfilter.hpp>
#include <dust2/unfilter.hpp>
#include <dust2/discrete/system.hpp>

namespace dust2 {
namespace r {

// TODO: this name must be changed!
template <typename T>
cpp11::sexp dust2_discrete_unfilter_alloc(cpp11::list r_pars,
                                          cpp11::sexp r_time_start,
                                          cpp11::sexp r_time,
                                          cpp11::sexp r_dt,
                                          cpp11::list r_data,
                                          cpp11::sexp r_n_particles,
                                          cpp11::sexp r_n_groups,
                                          cpp11::sexp r_index) {
  using rng_state_type = typename T::rng_state_type;

  auto n_particles = to_size(r_n_particles, "n_particles");
  auto n_groups = to_size(r_n_groups, "n_groups");
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
  auto seed = mcstate::random::r::as_rng_seed<rng_state_type>(R_NilValue);
  const auto deterministic = true;

  // Then allocate the system; this pulls together almost all the data
  // we need.  At this point we could have constructed the system out
  // of one that exists already on the R side, but I think that's
  // going to feel weirder overall.
  const auto system = dust2::dust_discrete<T>(shared, internal, time_start, dt, n_particles,
                                             seed, deterministic);
  const auto index = check_index(r_index, system.n_state(), "index");

  auto obj = new unfilter<dust_discrete<T>>(system, time_start, time, data, index);
  cpp11::external_pointer<unfilter<dust_discrete<T>>> ptr(obj, true, false);

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
