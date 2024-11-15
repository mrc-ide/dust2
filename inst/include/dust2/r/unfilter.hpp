#pragma once

#include <dust2/unfilter.hpp>
#include <monty/r/random.hpp>
#include <dust2/r/helpers.hpp>

namespace dust2 {
namespace r {

template <typename T>
cpp11::sexp dust2_unfilter_update_pars(cpp11::sexp ptr,
                                       cpp11::list r_pars,
                                       cpp11::sexp r_index_group) {
  auto *obj = dust2::r::safely_read_externalptr<unfilter<T>>(ptr, "unfilter_update_pars");
  const auto index_group = r_index_group == R_NilValue ? obj->sys.all_groups() :
    check_index(r_index_group, obj->sys.n_groups(), "index_group");
  update_pars(obj->sys, r_pars, index_group);
  return R_NilValue;
}

template <typename T>
cpp11::sexp dust2_unfilter_run(cpp11::sexp ptr, cpp11::sexp r_initial,
                               bool save_trajectories, bool adjoint,
                               cpp11::sexp r_index_state,
                               cpp11::sexp r_index_group,
                               bool preserve_particle_dimension,
                               bool preserve_group_dimension) {
  auto *obj = dust2::r::safely_read_externalptr<unfilter<T>>(ptr, "unfilter_run");
  const auto index_state = check_index(r_index_state, obj->sys.n_state(),
                                       "index_state");
  const auto index_group = r_index_group == R_NilValue ? obj->sys.all_groups() :
    check_index(r_index_group, obj->sys.n_groups(), "index_group");

  if (r_initial != R_NilValue) {
    set_state(obj->sys, cpp11::as_cpp<cpp11::list>(r_initial));
  }
  if (adjoint) {
    obj->run_adjoint(r_initial == R_NilValue, save_trajectories, index_state,
                     index_group);
  } else {
    obj->run(r_initial == R_NilValue, save_trajectories, index_state, index_group);
  }

  const auto n_groups = index_group.size();
  const auto n_particles = obj->sys.n_particles();
  const auto& ll = obj->last_log_likelihood();
  cpp11::writable::doubles ret(n_groups * n_particles);
  auto iter = REAL(ret);
  for (auto i : index_group) {
    iter = std::copy_n(ll.begin() + i * n_particles, n_particles, iter);
  }
  if (preserve_group_dimension && preserve_particle_dimension) {
    set_array_dims(ret, {n_particles, n_groups});
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_unfilter_last_trajectories(cpp11::sexp ptr,
                                             bool select_random_particle,
                                             bool preserve_particle_dimension,
                                             bool preserve_group_dimension) {
  auto *obj = dust2::r::safely_read_externalptr<unfilter<T>>(ptr, "unfilter_last_trajectories");

  const auto& trajectories = obj->last_trajectories();
  const auto& is_current = obj->last_trajectories_are_current();
  if (!tools::any(is_current)) {
    cpp11::stop("Trajectories are not current");
  }

  constexpr bool reorder = false; // never needed

  const auto n_state = trajectories.n_state(); // might be filtered
  const auto n_particles = trajectories.n_particles();
  const auto n_groups = trajectories.n_groups();
  const auto n_times = trajectories.n_times();

  const auto len = n_state * n_particles * n_groups * n_times;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  trajectories.export_state(REAL(ret), reorder, {});
  if (preserve_group_dimension && preserve_particle_dimension) {
    set_array_dims(ret, {n_state, n_particles, n_groups, n_times});
  } else if (preserve_group_dimension || preserve_particle_dimension) {
    set_array_dims(ret, {n_state, n_particles * n_groups, n_times});
  } else {
    set_array_dims(ret, {n_state * n_particles * n_groups, n_times});
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_unfilter_last_state(cpp11::sexp ptr,
                                      bool select_random_particle,
                                      bool preserve_particle_dimension,
                                      bool preserve_group_dimension) {
  auto *obj = dust2::r::safely_read_externalptr<unfilter<T>>(ptr, "unfilter_last_state");

  const auto& state = obj->sys.state();
  const auto n_state = obj->sys.n_state();
  const auto n_particles = obj->sys.n_particles();
  const auto n_groups = obj->last_index_group().size();

  const auto len = n_state * n_particles * n_groups;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  std::copy_n(state.begin(), len, REAL(ret));
  if (preserve_group_dimension && preserve_particle_dimension) {
    set_array_dims(ret, {n_state, n_particles, n_groups});
  } else if (preserve_group_dimension || preserve_particle_dimension) {
    set_array_dims(ret, {n_state, n_particles * n_groups});
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_discrete_unfilter_last_gradient(cpp11::sexp ptr,
                                                  bool preserve_particle_dimension,
                                                  bool preserve_group_dimension) {
  auto *obj = dust2::r::safely_read_externalptr<unfilter<T>>(ptr, "unfilter_last_gradient");
  const auto& is_current = obj->adjoint_is_current();
  if (!tools::any(is_current)) {
    cpp11::stop("System was not run with 'adjoint = TRUE'");
  }

  const auto n_gradient = obj->sys.packing_gradient().size();
  const auto n_particles = obj->sys.n_particles();
  const auto n_groups = obj->last_index_group().size();
  const auto len = n_gradient * n_particles * n_groups;

  if (n_particles > 1) {
    cpp11::stop("n_particles > 1 not yet supported; see mrc-5406");
  }
  cpp11::sexp ret = cpp11::writable::doubles(len);
  obj->last_gradient(REAL(ret));
  if (preserve_group_dimension && preserve_particle_dimension) {
    set_array_dims(ret, {n_gradient, n_particles, n_groups});
  } else if (preserve_group_dimension || preserve_particle_dimension) {
    set_array_dims(ret, {n_gradient, n_particles * n_groups});
  }
  return ret;
}

}
}
