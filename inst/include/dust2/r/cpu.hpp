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
                     cpp11::sexp r_seed,
                     cpp11::sexp r_deterministic) {
  using rng_state_type = typename T::rng_state_type;

  // Need to consider here the way that we are spreading out the
  // multiple parameter case _early_. For now we ignore that entirely
  // and assume we are doing just one.
  auto shared = T::build_shared(r_pars);
  auto internal = T::build_internal(r_pars);

  auto time = to_double(r_time, "time");
  auto dt = to_double(r_dt, "r_dt");
  auto n_particles = to_size(r_n_particles, "n_particles");
  auto seed = mcstate::random::r::as_rng_seed<rng_state_type>(r_seed);
  auto deterministic = to_bool(r_deterministic, "deterministic");

  auto obj = new dust_cpu<T>(shared, internal, time, dt, n_particles,
                             seed, deterministic);
  cpp11::external_pointer<dust_cpu<T>> ptr(obj, true, false);

  return cpp11::writable::list({ptr});
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
