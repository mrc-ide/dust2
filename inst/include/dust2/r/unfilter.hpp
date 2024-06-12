#pragma once

#include <mcstate/r/random.hpp>
#include <dust2/r/helpers.hpp>
#include <dust2/unfilter.hpp>

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

  auto obj = new unfilter<T>(system, time_start, time, data, index);
  cpp11::external_pointer<unfilter<T>> ptr(obj, true, false);

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

template <typename T>
cpp11::sexp dust2_discrete_unfilter_update_pars(cpp11::sexp ptr,
                                                cpp11::list r_pars,
                                                bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter<T>>>(ptr).get();
  update_pars(obj->sys, cpp11::as_cpp<cpp11::list>(r_pars), grouped);
  return R_NilValue;
}

template <typename T>
cpp11::sexp dust2_discrete_unfilter_run(cpp11::sexp ptr, cpp11::sexp r_initial,
                                        bool save_history, bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter<T>>>(ptr).get();
  if (r_initial != R_NilValue) {
    set_state(obj->sys, r_initial, grouped);
  }
  obj->run(r_initial == R_NilValue, save_history);

  const auto n_groups = obj->sys.n_groups();
  const auto n_particles = obj->sys.n_particles();
  cpp11::writable::doubles ret(n_groups * n_particles);
  obj->last_log_likelihood(REAL(ret));
  if (grouped && n_particles > 1) {
    set_array_dims(ret, {n_particles, n_groups});
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_discrete_unfilter_last_history(cpp11::sexp ptr, bool grouped) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter<T>>>(ptr).get();
  if (!obj->last_history_is_current()) {
    cpp11::stop("History is not current");
  }

  constexpr bool reorder = false; // never needed

  const auto& history = obj->last_history();
  const auto& dims = history.dims();
  // Could use destructured bind here in recent C++?
  const auto n_state = dims[0];
  const auto n_particles = dims[1];
  const auto n_groups = dims[2];
  const auto n_times = dims[3];
  const auto len = n_state * n_particles * n_groups * n_times;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  history.export_state(REAL(ret), reorder);
  if (grouped) {
    set_array_dims(ret, {n_state, n_particles, n_groups, n_times});
  } else {
    set_array_dims(ret, {n_state, n_particles * n_groups, n_times});
  }
  return ret;
}

}
}
