#pragma once

#include <dust2/r/helpers.hpp>
#include <mcstate/r/random.hpp>

#include <dust2/discrete.hpp>
#include <dust2/history.hpp>

namespace dust2 {
namespace r {

template <typename T>
SEXP dust2_discrete_alloc(cpp11::list r_pars,
                          cpp11::sexp r_time,
                          cpp11::sexp r_dt,
                          cpp11::sexp r_n_particles,
                          cpp11::sexp r_n_groups,
                          cpp11::sexp r_seed,
                          cpp11::sexp r_deterministic) {
  using rng_state_type = typename T::rng_state_type;

  const auto time = check_time(r_time, "time");
  const auto dt = check_dt(r_dt);

  const auto n_particles = to_size(r_n_particles, "n_particles");
  const auto n_groups = to_size(r_n_groups, "n_groups");

  const auto shared = build_shared<T>(r_pars, n_groups);
  // Later, we need one of these per thread
  const auto internal = build_internal<T>(shared);

  auto seed = mcstate::random::r::as_rng_seed<rng_state_type>(r_seed);
  auto deterministic = to_bool(r_deterministic, "deterministic");

  auto obj = new dust_discrete<T>(shared, internal, time, dt, n_particles,
                                  seed, deterministic);
  cpp11::external_pointer<dust_discrete<T>> ptr(obj, true, false);

  // Later, we'll export a bit more back from the system (in particular
  // systems need to provide information about how they organise
  // variables, ode systems report computed control, etc.
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

template <typename T>
SEXP dust2_discrete_run_steps(cpp11::sexp ptr, cpp11::sexp r_n_steps) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  auto n_steps = to_size(r_n_steps, "n_steps");
  obj->run_steps(n_steps);
  return R_NilValue;
}

template <typename T>
SEXP dust2_discrete_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  const auto time = check_time(r_time, "time");
  const auto curr = obj->time();
  if (time < curr) {
    cpp11::stop("Can't run to time %f, system already at time %f",
                time, curr);
  }
  obj->run_to_time(time);
  return R_NilValue;
}

template <typename T>
SEXP dust2_discrete_state(cpp11::sexp ptr, bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  cpp11::sexp ret = R_NilValue;
  const auto iter = obj->state().begin();
  if (grouped) {
    ret = export_array_n(iter,
                         {obj->n_state(), obj->n_particles(), obj->n_groups()});
  } else {
    ret = export_array_n(iter,
                         {obj->n_state(), obj->n_particles() * obj->n_groups()});
  }
  return ret;
}

template <typename T>
SEXP dust2_discrete_time(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  return cpp11::as_sexp(obj->time());
}

// If this is a grouped system then we will return a matrix with
// dimensions )(state x particle x group) and if we are ungrouped
// (state x particle); the difference is really only apparent in the
// case where we have a single group to drop.  In dust1 we had the
// option for (state x group) for simulations too.
//
// For now perhaps lets just ignore this detail and always return the
// 3d version as the code to do grouped/ungrouped switching is deep
// within one of the many in-progress branches and I can't find it
// yet.
//
// In the case where we have time, that goes in the last position
template <typename T>
SEXP dust2_discrete_set_state_initial(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  obj->set_state_initial();
  return R_NilValue;
}

template <typename T>
SEXP dust2_discrete_set_state(cpp11::sexp ptr, cpp11::sexp r_state, bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  set_state(*obj, r_state, grouped);
  return R_NilValue;
}

// Not really intended to be called by users, this just helps us test
// bookkeeping really.
template <typename T>
SEXP dust2_discrete_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  const auto n_particles = obj->n_particles();
  const auto n_groups = obj->n_groups();
  const int len = n_particles * n_groups;
  // We really should expect a matrix perhaps, but this test will
  // catch that confusingly at least.  Users won't actually call this.
  if (r_index.size() != len) {
    cpp11::stop("Expected an index of length %d", len);
  }
  std::vector<size_t> index;
  index.reserve(len);
  for (auto i : r_index) {
    if (i < 1 || i > len) {
      cpp11::stop("Expected 'index' values to lie in [1, %d]", n_particles);
    }
    index.push_back(i - 1);
  }
  obj->reorder(index.begin());
  return R_NilValue;
}

template <typename T>
SEXP dust2_discrete_rng_state(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  return rng_state_as_raw(obj->rng_state());
}

template <typename T>
SEXP dust2_discrete_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  using rng_state_type = typename T::rng_state_type;
  const auto n_streams = obj->n_particles();
  const auto rng_state =
    check_rng_state<rng_state_type>(r_rng_state, n_streams, "rng_state");
  obj->set_rng_state(rng_state);
  return R_NilValue;
}

template <typename T>
SEXP dust2_discrete_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  const auto time = check_time(r_time, "time");
  obj->set_time(time);
  return R_NilValue;
}

template <typename T>
SEXP dust2_discrete_update_pars(cpp11::sexp ptr, cpp11::list r_pars,
                           bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  update_pars(*obj, r_pars, grouped);
  return R_NilValue;
}

// This one exists to help push around the comparison part of things;
// it's not expected to be called often by users.
template <typename T>
SEXP dust2_discrete_compare_data(cpp11::sexp ptr,
                            cpp11::sexp r_data,
                            bool grouped) {
  using data_type = typename T::data_type;
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  const auto n_groups = obj->n_groups();
  std::vector<data_type> data;
  auto r_data_list = cpp11::as_cpp<cpp11::list>(r_data);
  if (grouped) {
    check_length(r_data_list, n_groups, "data");
    for (size_t i = 0; i < n_groups; ++i) {
      auto r_data_list_i = cpp11::as_cpp<cpp11::list>(r_data_list[i]);
      data.push_back(T::build_data(r_data_list_i));
    }
  } else {
    data.push_back(T::build_data(r_data_list));
  }

  cpp11::writable::doubles ret(obj->n_particles() * obj->n_groups());
  obj->compare_data(data.begin(), REAL(ret));
  if (grouped) {
    set_array_dims(ret, {obj->n_particles(), obj->n_groups()});
  }
  return ret;
}

template <typename T>
SEXP dust2_discrete_simulate(cpp11::sexp ptr,
                        cpp11::sexp r_times,
                        cpp11::sexp r_index,
                        bool grouped) {
  using real_type = typename T::real_type;
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_discrete<T>>>(ptr).get();
  const auto n_state = obj->n_state();
  const auto times = check_time_sequence(obj->time(), r_times, false, "time");
  const auto index = check_index(r_index, obj->n_state(), "index");
  const auto use_index = index.size() > 0;
  const auto n_state_save = use_index ? index.size() : n_state;
  const auto n_particles = obj->n_particles();
  const auto n_groups = obj->n_groups();
  const auto n_times = times.size();

  dust2::history<real_type> h(n_state_save, n_particles, n_groups, n_times);
  for (size_t i = 0; i < n_times; ++i) {
    obj->run_to_time(times[i]);
    if (use_index) {
      h.add_with_index(times[i], obj->state().begin(), index.begin(), n_state);
    } else {
      h.add(times[i], obj->state().begin());
    }
  }

  // There is an extra copy here vs using the R memory to back the
  // history.  That's an optimisation that would be fairly easy to
  // make later.
  const auto len = n_state_save * n_particles * n_groups * n_times;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  h.export_state(REAL(ret), false);
  if (grouped) {
    set_array_dims(ret, {n_state_save, n_particles, n_groups, n_times});
  } else {
    set_array_dims(ret, {n_state_save, n_particles * n_groups, n_times});
  }
  return ret;
}

}
}
