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

  // Then allocate the model; this pulls together almost all the data
  // we need.  At this point we could have constructed the model out
  // of one that exists already on the R side, but I think that's
  // going to feel weirder overall.
  const auto model = dust2::dust_cpu<T>(shared, internal, time_start, dt, n_particles,
                                 seed, deterministic);
  const auto index = check_index(r_index, model.n_state(), "index");

  auto obj = new unfilter<T>(model, time_start, time, data, index);
  cpp11::external_pointer<unfilter<T>> ptr(obj, true, false);

  cpp11::sexp r_n_state = cpp11::as_sexp(obj->model.n_state());
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

template <typename T>
cpp11::sexp dust2_cpu_unfilter_run(cpp11::sexp ptr, cpp11::sexp r_pars,
                                   cpp11::sexp r_initial, bool save_history,
                                   bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter<T>>>(ptr).get();
  if (r_pars != R_NilValue) {
    update_pars(obj->model, cpp11::as_cpp<cpp11::list>(r_pars), grouped);
  }
  if (r_initial != R_NilValue) {
    set_state(obj->model, r_initial, grouped);
  }
  obj->run(r_initial == R_NilValue, save_history);

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
cpp11::sexp dust2_cpu_unfilter_last_history(cpp11::sexp ptr, bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter<T>>>(ptr).get();
  if (!obj->last_history_is_current()) {
    cpp11::stop("History is not current");
  }

  const auto& dims = obj->last_history_dims();
  // Could use destructured bind here in recent C++?
  const auto n_state = dims[0];
  const auto n_particles = dims[1];
  const auto n_groups = dims[2];
  const auto n_times = dims[3];
  const auto len = n_state * n_particles * n_groups * n_times;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  obj->last_history(REAL(ret));
  if (grouped) {
    set_array_dims(ret, {n_state, n_particles, n_groups, n_times});
  } else {
    set_array_dims(ret, {n_state, n_particles * n_groups, n_times});
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
  // model, in a reasonable way.  We need to make this slightly easier
  // to do from mcstate really.  Expand the state to give all the
  // state required by the filter (n_groups streams worth) and the
  // model (n_groups * n_particles worth, though the last bit of
  // expansion could be done by the model itself instead?)
  //
  // There are two ways of sorting out the state here:
  //
  // 1. we could take the first n_groups states for the filter and the
  // remaining for the models.  This has the nice property that we can
  // expand the model state later if we support growing models
  // (mrc-5355).  However, it has the undesirable consequence that a
  // filter with multiple groups will stream differently to a filter
  // containing a the first group only.
  //
  // 2. we take each block of (1+n_particles) states for each group,
  // giving the first to the filter and the rest to the model.  This
  // means that we can change the number of groups without affecting
  // the results, though we can't change the number of particles as
  // easily.
  const auto n_groups_effective = grouped ? n_groups : 1;
  const auto n_streams = n_groups_effective * (n_particles + 1);
  const auto rng_state = mcstate::random::prng<rng_state_type>(n_streams, seed, deterministic).export_state();
  const auto rng_len = rng_state_type::size();
  rng_seed_type seed_filter;
  rng_seed_type seed_model;
  for (size_t i = 0; i < n_groups_effective; ++i) {
    const auto it = rng_state.begin() + i * rng_len * (n_particles + 1);
    seed_filter.insert(seed_filter.end(),
                       it, it + rng_len);
    seed_model.insert(seed_model.end(),
                      it + rng_len, it + rng_len * (n_particles + 1));
  }

  const auto model = dust2::dust_cpu<T>(shared, internal, time_start, dt, n_particles,
                                        seed_model, deterministic);

  const auto index = check_index(r_index, model.n_state(), "index");

  auto obj = new filter<T>(model, time_start, time, data, index, seed_filter);
  cpp11::external_pointer<filter<T>> ptr(obj, true, false);

  cpp11::sexp r_n_state = cpp11::as_sexp(obj->model.n_state());
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

template <typename T>
cpp11::sexp dust2_cpu_filter_run(cpp11::sexp ptr, cpp11::sexp r_pars,
                                 cpp11::sexp r_initial, bool save_history,
                                 bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  if (r_pars != R_NilValue) {
    update_pars(obj->model, cpp11::as_cpp<cpp11::list>(r_pars), grouped);
  }
  if (r_initial != R_NilValue) {
    set_state(obj->model, r_initial, grouped);
  }
  obj->run(r_initial == R_NilValue, save_history);

  cpp11::writable::doubles ret(obj->model.n_groups());
  obj->last_log_likelihood(REAL(ret));
  return ret;
}

// Can collapse with above
template <typename T>
cpp11::sexp dust2_cpu_filter_last_history(cpp11::sexp ptr, bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  if (!obj->last_history_is_current()) {
    cpp11::stop("History is not current");
  }

  const auto& dims = obj->last_history_dims();
  // Could use destructured bind here in recent C++?
  const auto n_state = dims[0];
  const auto n_particles = dims[1];
  const auto n_groups = dims[2];
  const auto n_times = dims[3];
  const auto len = n_state * n_particles * n_groups * n_times;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  obj->last_history(REAL(ret));
  if (grouped) {
    set_array_dims(ret, {n_state, n_particles, n_groups, n_times});
  } else {
    set_array_dims(ret, {n_state, n_particles * n_groups, n_times});
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_cpu_filter_rng_state(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  using rng_state_type = typename T::rng_state_type;

  // Undo the construction as above so that the rng state comes out in
  // the same format it goes in, as a single raw vector.
  const auto& state_filter = obj->rng_state();
  const auto& state_model = obj->model.rng_state();
  const auto n_particles = obj->model.n_particles();
  const auto n_groups = obj->model.n_groups();
  const auto n_state = rng_state_type::size();
  const auto n_bytes = sizeof(typename rng_state_type::int_type);
  const auto n_bytes_state = n_bytes * n_state;
  cpp11::writable::raws ret(n_bytes * (state_filter.size() + state_model.size()));
  for (size_t i = 0; i < n_groups; ++i) {
    std::memcpy(RAW(ret) + i * n_bytes_state * (n_particles + 1),
                state_filter.data() + i * n_state,
                n_bytes_state);
    std::memcpy(RAW(ret) + i * n_bytes_state * (n_particles + 1) + n_bytes_state,
                state_model.data() + i * n_state * n_particles,
                n_bytes_state * n_particles);
  }

  return ret;
}

}
}
