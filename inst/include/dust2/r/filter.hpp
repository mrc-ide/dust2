#pragma once

#include <dust2/filter.hpp>
#include <monty/r/random.hpp>
#include <dust2/r/helpers.hpp>

namespace dust2 {
namespace r {

template <typename T>
cpp11::sexp dust2_filter_update_pars(cpp11::sexp ptr,
                                     cpp11::list r_pars,
                                     cpp11::sexp r_index_group) {
  auto *obj = dust2::r::safely_read_externalptr<filter<T>>(ptr, "filter_update_pars");
  const auto index_group = r_index_group == R_NilValue ? obj->sys.all_groups() :
    check_index(r_index_group, obj->sys.n_groups(), "index_group");
  update_pars(obj->sys, r_pars, index_group);
  return R_NilValue;
}

template <typename T>
cpp11::sexp dust2_filter_run(cpp11::sexp ptr,
                             cpp11::sexp r_initial,
                             bool save_trajectories,
                             cpp11::sexp r_save_restart,
                             bool adjoint,
                             cpp11::sexp r_index_state,
                             cpp11::sexp r_index_group,
                             bool preserve_particle_dimension,
                             bool preserve_group_dimension) {
  auto *obj = dust2::r::safely_read_externalptr<filter<T>>(ptr, "filter_run");
  const auto index_state = check_index(r_index_state, obj->sys.n_state(),
                                       "index_state");
  const auto index_group = r_index_group == R_NilValue ? obj->sys.all_groups() :
    check_index(r_index_group, obj->sys.n_groups(), "index_group");
  using real_type = typename T::real_type;
  const auto save_restart = check_save_restart<real_type>(r_save_restart);

  if (r_initial != R_NilValue) {
    set_state(obj->sys, cpp11::as_cpp<cpp11::list>(r_initial));
  }
  obj->run(r_initial == R_NilValue, save_trajectories, save_restart, index_state, index_group);

  const auto& ll = obj->last_log_likelihood();
  cpp11::writable::doubles ret(index_group.size());
  auto iter = REAL(ret);
  for (auto i : index_group) {
    *iter = ll[i];
    iter++;
  }
  return ret;
}

// We might accept index_state here too later, allowing subsetting and
// validation of the index used.
template <typename T>
cpp11::sexp dust2_filter_last_trajectories(cpp11::sexp ptr,
                                           bool select_random_particle,
                                           bool preserve_particle_dimension,
                                           bool preserve_group_dimension) {
  auto *obj = dust2::r::safely_read_externalptr<filter<T>>(ptr, "filter_last_trajectories");

  const auto& trajectories = obj->last_trajectories();
  const auto& is_current = obj->last_trajectories_are_current();
  if (!tools::any(is_current)) {
    cpp11::stop("Trajectories are not current");
  }

  // We might relax this later, but will require some tools to work
  // with the output, really.
  constexpr bool reorder = true;

  const auto n_state = trajectories.n_state(); // might be filtered
  const auto n_particles = select_random_particle ? 1 : trajectories.n_particles();
  const auto n_groups = trajectories.n_groups();
  const auto n_times = trajectories.n_times();

  const auto& index_particle = select_random_particle ?
    obj->select_random_particle(trajectories.index_group()) : std::vector<size_t>{};

  const auto len = n_state * n_particles * n_groups * n_times;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  trajectories.export_state(REAL(ret), reorder, index_particle);
  if (select_random_particle) {
    if (preserve_group_dimension) {
      set_array_dims(ret, {n_state, n_particles * n_groups, n_times});
    } else {
      set_array_dims(ret, {n_state, n_particles * n_groups * n_times});
    }
  } else {
    if (preserve_group_dimension) {
      set_array_dims(ret, {n_state, n_particles, n_groups, n_times});
    } else {
      set_array_dims(ret, {n_state, n_particles * n_groups, n_times});
    }
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_filter_last_state(cpp11::sexp ptr,
                                    bool select_random_particle,
                                    bool preserve_particle_dimension,
                                    bool preserve_group_dimension) {
  auto *obj = dust2::r::safely_read_externalptr<filter<T>>(ptr, "filter_last_state");

  const auto& index_group = obj->last_index_group();

  const auto& state = obj->sys.state();
  const auto n_state = obj->sys.n_state();
  const auto n_particles = select_random_particle ? 1 : obj->sys.n_particles();
  const auto n_groups = index_group.size();

  const auto len = n_state * n_particles * n_groups;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  auto iter_dst = REAL(ret);

  if (select_random_particle) {
    const std::vector<size_t>& index_particle =
      obj->select_random_particle(index_group);
    const auto n_stride = n_state * obj->sys.n_particles();
    for (auto i : index_group) {
      const auto offset = i * n_stride + index_particle[i] * n_state;
      iter_dst = std::copy_n(state.begin() + offset, n_state, iter_dst);
    }
  } else {
    std::copy_n(state.begin(), len, iter_dst);
  }

  preserve_particle_dimension =
    preserve_particle_dimension && !select_random_particle;
  if (preserve_group_dimension && preserve_particle_dimension) {
    set_array_dims(ret, {n_state, n_particles, n_groups});
  } else if (preserve_group_dimension || preserve_particle_dimension) {
    set_array_dims(ret, {n_state, n_particles * n_groups});
  }
  return ret;
}

template <typename T>
cpp11::sexp dust2_filter_rng_state(cpp11::sexp ptr) {
  auto *obj = dust2::r::safely_read_externalptr<filter<T>>(ptr, "filter_rng_state");
  using rng_state_type = typename T::rng_state_type;

  // Undo the construction as above so that the rng state comes out in
  // the same format it goes in, as a single raw vector.
  const auto& state_filter = obj->rng_state();
  const auto& state_system = obj->sys.rng_state();
  const auto n_particles = obj->sys.n_particles();
  const auto n_groups = obj->sys.n_groups();
  const auto n_state = rng_state_type::size();
  const auto n_bytes = sizeof(typename rng_state_type::int_type);
  const auto n_bytes_state = n_bytes * n_state;
  cpp11::writable::raws ret(n_bytes * (state_filter.size() + state_system.size()));
  for (size_t i = 0; i < n_groups; ++i) {
    std::memcpy(RAW(ret) + i * n_bytes_state * (n_particles + 2),
                state_filter.data() + i * n_state * 2,
                n_bytes_state * 2);
    std::memcpy(RAW(ret) + i * n_bytes_state * (n_particles + 2) + 2 * n_bytes_state,
                state_system.data() + i * n_state * n_particles,
                n_bytes_state * n_particles);
  }

  return ret;
}

template <typename T>
cpp11::sexp dust2_filter_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  auto *obj = dust2::r::safely_read_externalptr<filter<T>>(ptr, "filter_set_rng_state");
  using rng_state_type = typename T::rng_state_type;
  using rng_int_type = typename rng_state_type::int_type;

  const auto n_particles = obj->sys.n_particles();
  const auto n_groups = obj->sys.n_groups();
  const auto n_streams = n_groups * (2 + n_particles);
  const auto n_state = rng_state_type::size();
  const auto rng_state =
    check_rng_state<rng_state_type>(r_rng_state, n_streams, "rng_state");

  std::vector<rng_int_type> state_filter(n_groups * n_state * 2);
  std::vector<rng_int_type> state_system(n_groups * n_particles * n_state);
  for (size_t i = 0; i < n_groups; ++i) {
    const auto src = rng_state.begin() + i * (2 + n_particles) * n_state;
    std::copy_n(src,
                n_state * 2,
                state_filter.begin() + i * n_state * 2);
    std::copy_n(src + 2 * n_state,
                n_state * n_particles,
                state_system.begin() + i * n_state * n_particles);
  }

  obj->set_rng_state(state_filter);
  obj->sys.set_rng_state(state_system);

  return R_NilValue;
}

}
}
