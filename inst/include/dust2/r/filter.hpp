#pragma once

#include <dust2/r/helpers.hpp>
#include <dust2/filter.hpp>

namespace dust2 {
namespace r {

// TODO: this name must be changed!
template <typename T>
cpp11::sexp dust2_cpu_unfilter_alloc(cpp11::list r_pars,
                                     cpp11::sexp r_time_start,
                                     cpp11::sexp r_time,
                                     cpp11::sexp r_dt,
                                     cpp11::list r_data,
                                     cpp11::sexp r_n_particles,
                                     cpp11::sexp r_n_groups) {
  using real_type = typename T::real_type;
  using rng_state_type = typename T::rng_state_type;

  auto n_particles = to_size(r_n_particles, "n_particles");
  auto n_groups = to_size(r_n_groups, "n_groups");
  const auto grouped = n_groups > 0;
  const auto time_start = check_time(r_time_start, "time_start");
  const auto time = check_time_sequence<real_type>(time_start, r_time, "time");
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

  // Then allocate the model; this pulls together almost all the data
  // we need.  At this point we could have constructed the model out
  // of one that exists already on the R side, but I think that's
  // going to feel weirder overall.
  const auto model = dust2::dust_cpu<T>(shared, internal, time_start, dt, n_particles,
                                 seed, deterministic);

  auto obj = new unfilter<T>(model, time_start, time, data);
  cpp11::external_pointer<unfilter<T>> ptr(obj, true, false);

  cpp11::sexp r_n_state = cpp11::as_sexp(obj->model.n_state());
  cpp11::sexp r_group_names = R_NilValue;
  if (grouped) {
    r_group_names = r_pars.attr("names");
  }
  cpp11::sexp r_grouped = cpp11::as_sexp(grouped);

  return cpp11::writable::list{ptr, r_n_state, r_grouped, r_group_names};
}

template <typename T>
cpp11::sexp dust2_cpu_unfilter_run(cpp11::sexp ptr, cpp11::sexp r_pars,
                                   cpp11::sexp r_initial, bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter<T>>>(ptr).get();
  if (r_pars != R_NilValue) {
    update_pars(obj->model, cpp11::as_cpp<cpp11::list>(r_pars), grouped);
  }
  if (r_initial != R_NilValue) {
    set_state(obj->model, r_initial, grouped);
  }
  obj->run(r_initial == R_NilValue);

  const auto n_groups = obj->model.n_groups();
  const auto n_particles = obj->model.n_particles();
  cpp11::writable::doubles ret(n_groups * n_particles);
  obj->last_log_likelihood(REAL(ret));
  if (grouped && n_particles > 1) {
    set_array_dims(ret, {n_particles, n_groups});
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_cpu_filter_alloc(cpp11::list r_pars,
                                   cpp11::sexp r_time_start,
                                   cpp11::sexp r_time,
                                   cpp11::sexp r_dt,
                                   cpp11::list r_data,
                                   cpp11::sexp r_n_particles,
                                   cpp11::sexp r_n_groups,
                                   cpp11::sexp r_seed) {
  using real_type = typename T::real_type;
  using rng_state_type = typename T::rng_state_type;
  using rng_seed_type = std::vector<typename rng_state_type::int_type>;

  auto n_particles = to_size(r_n_particles, "n_particles");
  auto n_groups = to_size(r_n_groups, "n_groups");
  const auto grouped = n_groups > 0;
  const auto time_start = check_time(r_time_start, "time_start");
  const auto time = check_time_sequence<real_type>(time_start, r_time, "time");
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
  // model, in a reasonable way.  We need to make this slightly easier
  // to do from mcstate really.  Expand the state to give all the
  // state required by the filter (n_groups streams worth) and the
  // model (n_groups * n_particles worth, though the last bit of
  // expansion could be done by the model itself instead?)
  const auto n_streams = n_groups * (n_particles + 1);
  const auto rng_state = mcstate::random::prng<rng_state_type>(n_streams, seed, deterministic).export_state();
  const auto rng_len = rng_state_type::size();
  const rng_seed_type seed_filter(rng_state.begin(),
                                  rng_state.begin() + rng_len * n_groups);
  const rng_seed_type seed_model(rng_state.begin() + rng_len * n_groups,
                                 rng_state.end());

  const auto model = dust2::dust_cpu<T>(shared, internal, time_start, dt, n_particles,
                                        seed_model, deterministic);

  auto obj = new filter<T>(model, time_start, time, data, seed_filter);
  cpp11::external_pointer<filter<T>> ptr(obj, true, false);

  cpp11::sexp r_n_state = cpp11::as_sexp(obj->model.n_state());
  cpp11::sexp r_group_names = R_NilValue;
  if (grouped) {
    r_group_names = r_pars.attr("names");
  }
  cpp11::sexp r_grouped = cpp11::as_sexp(grouped);

  return cpp11::writable::list{ptr, r_n_state, r_grouped, r_group_names};
}

template <typename T>
cpp11::sexp dust2_cpu_filter_run(cpp11::sexp ptr, cpp11::sexp r_pars,
                                 bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  if (r_pars != R_NilValue) {
    update_pars(obj->model, cpp11::as_cpp<cpp11::list>(r_pars), grouped);
  }
  obj->run();

  cpp11::writable::doubles ret(obj->model.n_groups());
  obj->last_log_likelihood(REAL(ret));
  return ret;
}

}
}
