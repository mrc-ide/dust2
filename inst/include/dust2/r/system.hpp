#pragma once

#include <dust2/r/helpers.hpp>
#include <dust2/history.hpp>

namespace dust2 {
namespace r {

template <typename T>
SEXP dust2_system_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
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
SEXP dust2_system_state(cpp11::sexp ptr, bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
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
SEXP dust2_system_time(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
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
SEXP dust2_system_set_state_initial(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  obj->set_state_initial();
  return R_NilValue;
}

template <typename T>
SEXP dust2_system_set_state(cpp11::sexp ptr, cpp11::sexp r_state, bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  set_state(*obj, r_state, grouped);
  return R_NilValue;
}

// Not really intended to be called by users, this just helps us test
// bookkeeping really.
template <typename T>
SEXP dust2_system_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
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
SEXP dust2_system_rng_state(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  return rng_state_as_raw(obj->rng_state());
}

template <typename T>
SEXP dust2_system_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  using rng_state_type = typename T::rng_state_type;
  const auto n_streams = obj->n_particles();
  const auto rng_state =
    check_rng_state<rng_state_type>(r_rng_state, n_streams, "rng_state");
  obj->set_rng_state(rng_state);
  return R_NilValue;
}

template <typename T>
SEXP dust2_system_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const auto time = check_time(r_time, "time");
  obj->set_time(time);
  return R_NilValue;
}

template <typename T>
SEXP dust2_system_update_pars(cpp11::sexp ptr, cpp11::list r_pars,
                              bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  update_pars(*obj, r_pars, grouped);
  return R_NilValue;
}

// This one exists to help push around the comparison part of things;
// it's not expected to be called often by users.
template <typename T>
SEXP dust2_system_compare_data(cpp11::sexp ptr,
                               cpp11::sexp r_data,
                               bool grouped) {
  using system_type = typename T::system_type;
  using data_type = typename T::data_type;

  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  const auto n_groups = obj->n_groups();
  std::vector<data_type> data;
  auto r_data_list = cpp11::as_cpp<cpp11::list>(r_data);
  if (grouped) {
    check_length(r_data_list, n_groups, "data");
    for (size_t i = 0; i < n_groups; ++i) {
      auto r_data_list_i = cpp11::as_cpp<cpp11::list>(r_data_list[i]);
      data.push_back(system_type::build_data(r_data_list_i));
    }
  } else {
    data.push_back(system_type::build_data(r_data_list));
  }

  cpp11::writable::doubles ret(obj->n_particles() * obj->n_groups());
  obj->compare_data(data.begin(), REAL(ret));
  if (grouped) {
    set_array_dims(ret, {obj->n_particles(), obj->n_groups()});
  }
  return ret;
}

template <typename T>
SEXP dust2_system_simulate(cpp11::sexp ptr,
                           cpp11::sexp r_times,
                           cpp11::sexp r_index,
                           bool grouped) {
  using real_type = typename T::real_type;
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
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
