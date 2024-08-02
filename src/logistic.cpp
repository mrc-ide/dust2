// Generated by dust2 (version 0.1.0) - do not edit

#include <dust2/common.hpp>
#include <numeric>

// [[dust2::class(logistic)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(n)]]
// [[dust2::parameter(r)]]
// [[dust2::parameter(K)]]
class logistic {
public:
  logistic() = delete;

  using real_type = double;

  struct shared_state {
    size_t n;
    std::vector<real_type> r;
    std::vector<real_type> K;
  };

  using internal_state = dust2::no_internal_state;
  using data_type = dust2::no_data;
  using rng_state_type = mcstate::random::generator<real_type>;

  static size_t size_state(const shared_state& shared) {
    return shared.n;
  }

  static size_t size_output(const shared_state& shared) {
    return 1;
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state) {
    for (size_t i = 0; i < shared.n; ++i) {
      state[i] = 1;
    }
  }

  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
    for (size_t i = 0; i < shared.n; ++i) {
      state_deriv[i] = shared.r[i] * state[i] * (1 - state[i] / shared.K[i]);
    }
  }

  static void output(real_type time,
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     real_type * output) {
    // We will change this to use a delay (e.g., growth over last
    // period) to give the history a good workout right away.
    output[0] = std::accumulate(state, state + shared.n, 0);
  }

  static shared_state build_shared(cpp11::list pars) {
    const auto n = dust2::r::read_size(pars, "n");
    std::vector<real_type> r(n);
    std::vector<real_type> K(n);
    dust2::r::read_real_vector(pars, n, r.data(), "r", true);
    dust2::r::read_real_vector(pars, n, K.data(), "K", true);
    return shared_state{n, r, K};
  }

  static internal_state build_internal(const shared_state& shared) {
    return internal_state{};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
    dust2::r::read_real_vector(pars, shared.n, shared.r.data(), "r", false);
    dust2::r::read_real_vector(pars, shared.n, shared.K.data(), "K", false);
  }

  static void update_internal(const shared_state& shared,
                              internal_state& internal) {
  }

  static auto zero_every(const shared_state& shared) {
    return dust2::zero_every_type<real_type>();
  }
};

#include <cpp11.hpp>
#include <dust2/r/continuous/system.hpp>

[[cpp11::register]]
SEXP dust2_system_logistic_alloc(cpp11::list r_pars, cpp11::sexp r_time, cpp11::list r_ode_control, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_seed, cpp11::sexp r_deterministic, cpp11::sexp r_n_threads) {
  return dust2::r::dust2_continuous_alloc<logistic>(r_pars, r_time, r_ode_control, r_n_particles, r_n_groups, r_seed, r_deterministic, r_n_threads);
}

[[cpp11::register]]
SEXP dust2_system_logistic_internals(cpp11::sexp ptr, bool include_coefficients) {
  return dust2::r::dust2_system_internals<dust2::dust_continuous<logistic>>(ptr, include_coefficients);
}
[[cpp11::register]]
SEXP dust2_system_logistic_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_run_to_time<dust2::dust_continuous<logistic>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_logistic_state(cpp11::sexp ptr, cpp11::sexp r_index_state, cpp11::sexp r_index_particle, cpp11::sexp r_index_group, bool preserve_group_dimension) {
  return dust2::r::dust2_system_state<dust2::dust_continuous<logistic>>(ptr, r_index_state, r_index_particle, r_index_group, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_system_logistic_time(cpp11::sexp ptr) {
  return dust2::r::dust2_system_time<dust2::dust_continuous<logistic>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_logistic_set_state_initial(cpp11::sexp ptr) {
  return dust2::r::dust2_system_set_state_initial<dust2::dust_continuous<logistic>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_logistic_set_state(cpp11::sexp ptr, cpp11::sexp r_state) {
  return dust2::r::dust2_system_set_state<dust2::dust_continuous<logistic>>(ptr, r_state);
}

[[cpp11::register]]
SEXP dust2_system_logistic_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  return dust2::r::dust2_system_reorder<dust2::dust_continuous<logistic>>(ptr, r_index);
}

[[cpp11::register]]
SEXP dust2_system_logistic_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_system_rng_state<dust2::dust_continuous<logistic>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_logistic_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  return dust2::r::dust2_system_set_rng_state<dust2::dust_continuous<logistic>>(ptr, r_rng_state);
}

[[cpp11::register]]
SEXP dust2_system_logistic_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_set_time<dust2::dust_continuous<logistic>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_logistic_update_pars(cpp11::sexp ptr, cpp11::list pars) {
  return dust2::r::dust2_system_update_pars<dust2::dust_continuous<logistic>>(ptr, pars);
}

[[cpp11::register]]
SEXP dust2_system_logistic_simulate(cpp11::sexp ptr, cpp11::sexp r_times, cpp11::sexp r_index, bool preserve_group_dimension) {
  return dust2::r::dust2_system_simulate<dust2::dust_continuous<logistic>>(ptr, r_times, r_index, preserve_group_dimension);
}
