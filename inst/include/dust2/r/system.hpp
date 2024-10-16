#pragma once

#include <dust2/r/helpers.hpp>
#include <dust2/history.hpp>

namespace dust2 {
namespace r {

template <typename T>
void check_errors(T* obj, const char *action) {
  // It would be nice to throw a classed error here, but this is not
  // easy; see https://github.com/r-lib/cpp11/issues/250
  if (obj->errors_pending()) {
    cpp11::stop("Can't currently %s this system: errors are pending", action);
  }
}

template <typename T>
SEXP dust2_system_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  check_errors(obj, "run");
  const auto time = check_time(r_time, "time");
  const auto curr = obj->time();
  if (time < curr) {
    cpp11::stop("Can't run to time %f, system already at time %f",
                time, curr);
  }
  obj->run_to_time(time, obj->all_groups());
  return R_NilValue;
}

template <typename T>
SEXP dust2_system_state(cpp11::sexp ptr, cpp11::sexp r_index_state,
                        cpp11::sexp r_index_particle,
                        cpp11::sexp r_index_group,
                        bool preserve_particle_dimension,
                        bool preserve_group_dimension) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  check_errors(obj, "get state from");
  const auto n_state = obj->n_state();
  const auto n_particles = obj->n_particles();
  const auto n_groups = obj->n_groups();

  cpp11::sexp ret = R_NilValue;
  const auto iter_src = obj->state().begin();
  bool everything = r_index_state == R_NilValue &&
    r_index_particle == R_NilValue &&
    r_index_group == R_NilValue;
  if (everything) {
    if (preserve_group_dimension && preserve_particle_dimension) {
      ret = export_array_n(iter_src, {n_state, n_particles, n_groups});
    } else if (preserve_group_dimension || preserve_particle_dimension) {
      ret = export_array_n(iter_src, {n_state, n_particles * n_groups});
    } else {
      ret = export_array_n(iter_src, {n_state * n_particles * n_groups});
    }
  } else {
    // TODO: allow dropping these dimensions here, probably in another
    // ticket?
    if (!preserve_group_dimension && r_index_group != R_NilValue) {
      cpp11::stop("Can't provide 'index_group' for a non-grouped system");
    }
    const auto index_state =
      check_index(r_index_state, n_state, "index_state");
    const auto index_particle =
      check_index(r_index_particle, n_particles, "index_particle");
    const auto index_group =
      check_index(r_index_group, n_groups, "index_group");

    // This is surprisingly icky to do, but this will do a subsetting copy
    const auto n_state_save =
      index_state.empty() ? n_state : index_state.size();
    const auto n_particle_save =
      index_particle.empty() ? n_particles : index_particle.size();
    const auto n_group_save =
      index_group.empty() ? n_groups : index_group.size();
    cpp11::writable::doubles d(n_state_save * n_particle_save * n_group_save);
    double* it_dst = REAL(d);
    if (preserve_group_dimension && preserve_particle_dimension) {
      set_array_dims(d, {n_state_save, n_particle_save, n_group_save});
    } else if (preserve_group_dimension || preserve_particle_dimension) {
      set_array_dims(d, {n_state_save, n_particle_save * n_group_save});
    }

    for (size_t i = 0; i < n_group_save; ++i) {
      const auto ii = index_group.empty() ? i : index_group[i];
      auto it_i = iter_src + ii * n_state * n_particles;
      for (size_t j = 0; j < n_particle_save; ++j) {
        const auto jj = index_particle.empty() ? j : index_particle[j];
        const auto it_j = it_i + jj * n_state;
        if (index_state.empty()) {
          it_dst = std::copy_n(it_j, n_state, it_dst);
        } else {
          for (size_t k = 0; k < n_state_save; ++k, ++it_dst) {
            const auto kk = index_state[k];
            *it_dst = *(it_j + kk);
          }
        }
      }
    }
    ret = d;
  }
  return ret;
}

template <typename T>
SEXP dust2_system_time(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  check_errors(obj, "get time from");
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
  obj->set_state_initial(obj->all_groups());
  return R_NilValue;
}

template <typename T>
SEXP dust2_system_set_state(cpp11::sexp ptr,
                            cpp11::list r_state) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  set_state(*obj, r_state);
  return R_NilValue;
}

// Not really intended to be called by users, this just helps us test
// bookkeeping really.
template <typename T>
SEXP dust2_system_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  check_errors(obj, "reorder");
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
  obj->reorder(index.begin(), obj->all_groups());
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
  check_errors(obj, "set time for");
  const auto time = check_time(r_time, "time");
  obj->set_time(time);
  return R_NilValue;
}

template <typename T>
SEXP dust2_system_update_pars(cpp11::sexp ptr, cpp11::list r_pars) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  update_pars(*obj, r_pars, obj->all_groups());
  return R_NilValue;
}

// This one exists to help push around the comparison part of things;
// it's not expected to be called often by users.
template <typename T>
SEXP dust2_system_compare_data(cpp11::sexp ptr,
                               cpp11::list r_data,
                               bool preserve_particle_dimension,
                               bool preserve_group_dimension) {
  using system_type = typename T::system_type;
  using data_type = typename system_type::data_type;

  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  check_errors(obj, "compare data for");
  const auto n_groups = obj->n_groups();
  const auto& shared = obj->shared();

  std::vector<data_type> data;
  check_length(r_data, n_groups, "data");
  for (size_t i = 0; i < n_groups; ++i) {
    auto r_data_i = cpp11::as_cpp<cpp11::list>(r_data[i]);
    data.push_back(system_type::build_data(r_data_i, shared[i]));
  }

  cpp11::writable::doubles ret(obj->n_particles() * obj->n_groups());
  obj->compare_data(data.begin(), obj->all_groups(), REAL(ret));
  if (preserve_group_dimension && preserve_particle_dimension) {
    set_array_dims(ret, {obj->n_particles(), obj->n_groups()});
  }
  return ret;
}

template <typename T>
SEXP dust2_system_simulate(cpp11::sexp ptr,
                           cpp11::sexp r_times,
                           cpp11::sexp r_index_state,
                           bool preserve_particle_dimension,
                           bool preserve_group_dimension) {
  using real_type = typename T::real_type;
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<T>>(ptr).get();
  check_errors(obj, "simulate");
  const auto n_state = obj->n_state();
  const auto times = check_time_sequence(obj->time(), r_times, false, "time");
  const auto index_state =
    check_index(r_index_state, obj->n_state(), "index_state");
  const auto index_group = obj->all_groups();
  const auto use_index_state = index_state.size() > 0;
  const auto n_state_save = use_index_state ? index_state.size() : n_state;
  const auto n_particles = obj->n_particles();
  const auto n_groups = index_group.size();
  const auto n_times = times.size();

  dust2::history<real_type> h(n_state, n_particles, n_groups, n_times);
  h.set_index_and_reset(index_state, index_group);
  for (size_t i = 0; i < n_times; ++i) {
    obj->run_to_time(times[i], index_group);
    h.add(times[i], obj->state().begin());
  }

  // There is an extra copy here vs using the R memory to back the
  // history.  That's an optimisation that would be fairly easy to
  // make later.
  const auto len = n_state_save * n_particles * n_groups * n_times;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  const auto reorder = false;
  h.export_state(REAL(ret), reorder, {});
  if (preserve_group_dimension && preserve_particle_dimension) {
    set_array_dims(ret, {n_state_save, n_particles, n_groups, n_times});
  } else if (preserve_group_dimension || preserve_particle_dimension) {
    set_array_dims(ret, {n_state_save, n_particles * n_groups, n_times});
  } else {
    set_array_dims(ret, {n_state_save * n_particles * n_groups, n_times});
  }
  return ret;
}

}
}
