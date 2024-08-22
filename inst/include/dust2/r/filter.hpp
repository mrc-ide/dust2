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
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  const auto index_group = r_index_group == R_NilValue ? obj->sys.all_groups() :
    check_index(r_index_group, obj->sys.n_groups(), "index_group");
  update_pars(obj->sys, r_pars, index_group);
  return R_NilValue;
}

template <typename T>
cpp11::sexp dust2_filter_run(cpp11::sexp ptr, cpp11::sexp r_initial,
                             bool save_history,
                             cpp11::sexp r_index_group,
                             bool preserve_group_dimension) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  const auto index_group = r_index_group == R_NilValue ? obj->sys.all_groups() :
    check_index(r_index_group, obj->sys.n_groups(), "index_group");
  if (r_initial != R_NilValue) {
    set_state(obj->sys, r_initial, preserve_group_dimension, index_group);
  }
  obj->run(r_initial == R_NilValue, save_history, index_group);

  const auto& ll = obj->last_log_likelihood();
  cpp11::writable::doubles ret(index_group.size());
  auto iter = REAL(ret);
  for (auto i : index_group) {
    *iter = ll[i];
    iter++;
  }
  return ret;
}

// Can collapse with above
template <typename T>
cpp11::sexp dust2_filter_last_history(cpp11::sexp ptr,
                                      cpp11::sexp r_index_group,
                                      bool select_random_particle,
                                      bool preserve_group_dimension) {
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  const auto index_group = r_index_group == R_NilValue ? obj->sys.all_groups() :
    check_index(r_index_group, obj->sys.n_groups(), "index_group");
  const auto& is_current = obj->last_history_is_current();
  for (auto i : index_group) {
    if (!is_current[i]) {
      if (!tools::any(is_current)) {
        cpp11::stop("History is not current");
      } else {
        cpp11::stop("History for group '%d' is not current",
                    static_cast<int>(i + 1));
      }
    }
  }

  // We might relax this later, but will require some tools to work
  // with the output, really.
  constexpr bool reorder = true;

  const auto& history = obj->last_history();

  const auto n_state = history.n_state(); // might be filtered
  const auto n_particles = select_random_particle ? 1 : obj->sys.n_particles();
  const auto n_groups = index_group.size();
  const auto n_times = history.n_times();

  std::vector<size_t> index_particle;
  if (select_random_particle) {
    index_particle.resize(n_groups);
    for (auto i : index_group) {
      index_particle[i] = obj->select_random_particle(i);
    }
  }

  const auto len = n_state * n_particles * n_groups * n_times;
  cpp11::sexp ret = cpp11::writable::doubles(len);
  history.export_state(REAL(ret), reorder, index_group, index_particle);
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
cpp11::sexp dust2_filter_rng_state(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
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
    std::memcpy(RAW(ret) + i * n_bytes_state * (n_particles + 1),
                state_filter.data() + i * n_state,
                n_bytes_state);
    std::memcpy(RAW(ret) + i * n_bytes_state * (n_particles + 1) + n_bytes_state,
                state_system.data() + i * n_state * n_particles,
                n_bytes_state * n_particles);
  }

  return ret;
}

template <typename T>
cpp11::sexp dust2_filter_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<filter<T>>>(ptr).get();
  using rng_state_type = typename T::rng_state_type;
  using rng_int_type = typename rng_state_type::int_type;

  const auto n_particles = obj->sys.n_particles();
  const auto n_groups = obj->sys.n_groups();
  const auto n_streams = n_groups * (1 + n_particles);
  const auto n_state = rng_state_type::size();
  const auto rng_state =
    check_rng_state<rng_state_type>(r_rng_state, n_streams, "rng_state");

  std::vector<rng_int_type> state_filter(n_groups * n_state);
  std::vector<rng_int_type> state_system(n_groups * n_particles * n_state);
  for (size_t i = 0; i < n_groups; ++i) {
    const auto src = rng_state.begin() + i * (1 + n_particles) * n_state;
    std::copy_n(src,
                n_state,
                state_filter.begin() + i * n_state);
    std::copy_n(src + n_state,
                n_state * n_particles,
                state_system.begin() + i * n_state * n_particles);
  }

  obj->set_rng_state(state_filter);
  obj->sys.set_rng_state(state_system);

  return R_NilValue;
}

}
}
