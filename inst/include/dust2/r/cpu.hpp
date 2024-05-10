#pragma once

#include <dust2/r/helpers.hpp>
#include <mcstate/r/random.hpp>

#include <dust2/cpu.hpp>

namespace dust2 {
namespace r {

template <typename T>
SEXP dust2_cpu_alloc(cpp11::list r_pars,
                     cpp11::sexp r_time,
                     cpp11::sexp r_dt,
                     cpp11::sexp r_n_particles,
                     cpp11::sexp r_n_groups,
                     cpp11::sexp r_seed,
                     cpp11::sexp r_deterministic) {
  using shared_state = typename T::shared_state;
  using internal_state = typename T::internal_state;
  using rng_state_type = typename T::rng_state_type;

  const auto time = check_time(r_time);
  const auto dt = check_dt(r_dt);

  auto n_particles = to_size(r_n_particles, "n_particles");
  auto n_groups = to_size(r_n_groups, "n_groups");

  std::vector<shared_state> shared;
  std::vector<internal_state> internal;

  size_t size = 0;
  cpp11::sexp group_names = R_NilValue;
  const auto grouped = n_groups > 0;
  if (!grouped) {
    shared.push_back(T::build_shared(r_pars));
    internal.push_back(T::build_internal(r_pars));
    size = T::size(shared[0]);
  } else {
    if (r_pars.size() != static_cast<int>(n_groups)) {
      cpp11::stop("Expected 'pars' to have length %d to match 'n_groups'",
                  static_cast<int>(n_groups));
    }
    for (size_t i = 0; i < n_groups; ++i) {
      shared.push_back(T::build_shared(r_pars[i]));
      internal.push_back(T::build_internal(r_pars[i]));
      const auto size_i = T::size(shared[i]);
      if (i == 0) {
        size = size_i;
      } else if (size_i != size) {
        cpp11::stop("Expected state length for group %d to be %d, but it was %d",
                    i + 1, size, size_i);
      }
    }
    group_names = r_pars.attr("names");
  }

  auto seed = mcstate::random::r::as_rng_seed<rng_state_type>(r_seed);
  auto deterministic = to_bool(r_deterministic, "deterministic");

  auto obj = new dust_cpu<T>(shared, internal, time, dt, n_particles,
                             seed, deterministic);
  cpp11::external_pointer<dust_cpu<T>> ptr(obj, true, false);

  // Later, we'll export a bit more back from the model (in particular
  // models need to provide information about ~how they organise
  // variables~, ode models report computed control, etc.

  cpp11::sexp r_size = cpp11::as_sexp(size);
  cpp11::sexp r_grouped = cpp11::as_sexp(grouped);

  return cpp11::writable::list{ptr, r_size, r_grouped, group_names};
}

template <typename T>
SEXP dust2_cpu_run_steps(cpp11::sexp ptr, cpp11::sexp r_n_steps) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  auto n_steps = to_size(r_n_steps, "n_steps");
  obj->run_steps(n_steps);
  return R_NilValue;
}

template <typename T>
SEXP dust2_cpu_state(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  // return to_matrix(obj->state(), obj->n_state, obj->n_particles);
  return cpp11::as_sexp(obj->state());
}

template <typename T>
SEXP dust2_cpu_time(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  return cpp11::as_sexp(obj->time());
}

template <typename T>
SEXP dust2_cpu_set_state_initial(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  obj->set_state_initial();
  return R_NilValue;
}

template <typename T>
SEXP dust2_cpu_set_state(cpp11::sexp ptr, cpp11::sexp r_state) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  obj->set_state(REAL(r_state));
  return R_NilValue;
}

template <typename T>
SEXP dust2_cpu_rng_state(cpp11::sexp ptr) {
  using rng_state_type = typename T::rng_state_type;
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();

  const auto state = obj->rng_state();
  const auto len = sizeof(typename rng_state_type::int_type) * state.size();
  cpp11::writable::raws ret(len);
  std::memcpy(RAW(ret), state.data(), len);
  return ret;
}

template <typename T>
SEXP dust2_cpu_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  const auto time = check_time(r_time);
  obj->set_time(time);
  return R_NilValue;
}

template <typename T>
SEXP dust2_cpu_update_pars(cpp11::sexp ptr, cpp11::list r_pars,
                           bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  if (grouped) {
    const auto n_groups = obj->n_groups();
    if (r_pars.size() != static_cast<int>(n_groups)) {
      cpp11::stop("Expected 'pars' to have length %d to match 'n_groups'",
                  static_cast<int>(n_groups));
    }
    for (size_t i = 0; i < n_groups; ++i) {
      obj->update_shared(i, [&] (auto& shared) {
                              T::update_shared(r_pars[i], shared);
                            });
    }
  } else {
    obj->update_shared(0, [&] (auto& shared) {
                            T::update_shared(r_pars, shared);
                          });
  }
  return R_NilValue;
}

}
}
