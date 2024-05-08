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

  auto time = to_double(r_time, "time");
  auto dt = to_double(r_dt, "r_dt");
  auto n_particles = to_size(r_n_particles, "n_particles");
  auto n_groups = to_size(r_n_groups, "n_groups");

  std::vector<shared_state> shared;
  std::vector<internal_state> internal;

  size_t size = 0;
  cpp11::sexp group_names = R_NilValue;
  if (n_groups == 0) {
    shared.push_back(T::build_shared(r_pars));
    internal.push_back(T::build_internal(r_pars));
    size = T::size(shared[0]);
  } else {
    if (r_pars.size() != static_cast<int>(n_groups)) {
      cpp11::stop("Expected 'pars' to have length %d to match n_groups",
                  static_cast<int>(n_groups));
    }
    for (size_t i = 0; i < n_groups; ++i) {
      shared.push_back(T::build_shared(r_pars[i]));
      internal.push_back(T::build_internal(r_pars[i]));
      if (i == 0) {
        size = T::size(shared[i]);
      } else if (T::size(shared[i]) != size) {
        cpp11::stop("Invalid size generated for group %d", i + 1);
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
  // models need to provide information about how they organise
  // variables, ode models report computed control, etc.

  return cpp11::writable::list{ptr, cpp11::as_sexp(size), group_names};
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

}
}
