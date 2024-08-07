// Generated by dust2 (version 0.1.0) - do not edit

#include <dust2/common.hpp>

// [[dust2::class(sirode)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(I0)]]
// [[dust2::parameter(N)]]
// [[dust2::parameter(beta)]]
// [[dust2::parameter(gamma)]]
// [[dust2::parameter(exp_noise)]]
class sirode {
public:
  sirode() = delete;

  using real_type = double;

  struct shared_state {
    real_type N;
    real_type I0;
    real_type beta;
    real_type gamma;
    real_type exp_noise;
  };

  using internal_state = dust2::no_internal_state;

  struct data_type {
    real_type incidence;
  };

  using rng_state_type = mcstate::random::generator<real_type>;

  static size_t size_state(const shared_state& shared) {
    return 5;
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    state_next[0] = shared.N - shared.I0;
    state_next[1] = shared.I0;
    state_next[2] = 0;
    state_next[3] = 0;
    state_next[4] = 0;
  }

  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
    const auto S = state[0];
    const auto I = state[1];
    const auto rate_SI = shared.beta * S * I / shared.N;
    const auto rate_IR = shared.gamma * I;
    state_deriv[0] = -rate_SI;
    state_deriv[1] = rate_SI - rate_IR;
    state_deriv[2] = rate_IR;
    state_deriv[3] = rate_SI;
    state_deriv[4] = rate_SI;
  }

  static shared_state build_shared(cpp11::list pars) {
    const real_type I0 = dust2::r::read_real(pars, "I0", 10);
    const real_type N = dust2::r::read_real(pars, "N", 1000);
    const real_type beta = dust2::r::read_real(pars, "beta", 0.2);
    const real_type gamma = dust2::r::read_real(pars, "gamma", 0.1);
    const real_type exp_noise = dust2::r::read_real(pars, "exp_noise", 1e6);
    return shared_state{N, I0, beta, gamma, exp_noise};
  }

  static internal_state build_internal(const shared_state& shared) {
    return internal_state{};
  }

  // This is the bit that we'll use to do fast parameter updating, and
  // we'll guarantee somewhere that the size does not change.
  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.I0 = dust2::r::read_real(pars, "I0", shared.I0);
    shared.beta = dust2::r::read_real(pars, "beta", shared.beta);
    shared.gamma = dust2::r::read_real(pars, "gamma", shared.gamma);
  }

  // This is a reasonable default implementation in the no-internal
  // case
  static void update_internal(const shared_state& shared,
                              internal_state& internal) {
  }

  static auto zero_every(const shared_state& shared) {
    return dust2::zero_every_type<real_type>{{1, {4}}}; // zero[1] = {4};
  }
};

#include <cpp11.hpp>
#include <dust2/r/continuous/system.hpp>

[[cpp11::register]]
SEXP dust2_system_sirode_alloc(cpp11::list r_pars, cpp11::sexp r_time, cpp11::list r_ode_control, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_seed, cpp11::sexp r_deterministic, cpp11::sexp r_n_threads) {
  return dust2::r::dust2_continuous_alloc<sirode>(r_pars, r_time, r_ode_control, r_n_particles, r_n_groups, r_seed, r_deterministic, r_n_threads);
}

[[cpp11::register]]
SEXP dust2_system_sirode_internals(cpp11::sexp ptr, bool include_coefficients) {
  return dust2::r::dust2_system_internals<dust2::dust_continuous<sirode>>(ptr, include_coefficients);
}
[[cpp11::register]]
SEXP dust2_system_sirode_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_run_to_time<dust2::dust_continuous<sirode>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_sirode_state(cpp11::sexp ptr, cpp11::sexp r_index_state, cpp11::sexp r_index_particle, cpp11::sexp r_index_group, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_system_state<dust2::dust_continuous<sirode>>(ptr, r_index_state, r_index_particle, r_index_group, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_system_sirode_time(cpp11::sexp ptr) {
  return dust2::r::dust2_system_time<dust2::dust_continuous<sirode>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_sirode_set_state_initial(cpp11::sexp ptr) {
  return dust2::r::dust2_system_set_state_initial<dust2::dust_continuous<sirode>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_sirode_set_state(cpp11::sexp ptr, cpp11::sexp r_state, bool preserve_group_dimension) {
  return dust2::r::dust2_system_set_state<dust2::dust_continuous<sirode>>(ptr, r_state, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_system_sirode_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  return dust2::r::dust2_system_reorder<dust2::dust_continuous<sirode>>(ptr, r_index);
}

[[cpp11::register]]
SEXP dust2_system_sirode_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_system_rng_state<dust2::dust_continuous<sirode>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_sirode_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  return dust2::r::dust2_system_set_rng_state<dust2::dust_continuous<sirode>>(ptr, r_rng_state);
}

[[cpp11::register]]
SEXP dust2_system_sirode_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_set_time<dust2::dust_continuous<sirode>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_sirode_update_pars(cpp11::sexp ptr, cpp11::list pars) {
  return dust2::r::dust2_system_update_pars<dust2::dust_continuous<sirode>>(ptr, pars);
}

[[cpp11::register]]
SEXP dust2_system_sirode_simulate(cpp11::sexp ptr, cpp11::sexp r_times, cpp11::sexp r_index_state, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_system_simulate<dust2::dust_continuous<sirode>>(ptr, r_times, r_index_state, preserve_particle_dimension, preserve_group_dimension);
}
