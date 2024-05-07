#include <dust2/common.hpp>

namespace {
inline double with_default(double default_value, cpp11::sexp value) {
  return value == R_NilValue ? default_value : cpp11::as_cpp<double>(value);
}
}

class sir {
public:
  sir() = delete;

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

  static size_t size(const shared_state& shared) {
    return 5;
  }

  static void initial(real_type time,
                      real_type dt,
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

  // The main update function, converting state to state_next
  static void update(real_type time,
                     real_type dt,
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     rng_state_type& rng_state,
                     real_type * state_next) {
    const auto S = state[0];
    const auto I = state[1];
    const auto R = state[2];
    const auto cases_cumul = state[3];
    // const auto cases_inc = state[4];
    const auto p_SI = 1 - mcstate::math::exp(-shared.beta * I / shared.N * dt);
    const auto p_IR = 1 - mcstate::math::exp(-shared.gamma * dt);
    const auto n_SI = mcstate::random::binomial<real_type>(rng_state, S, p_SI);
    const auto n_IR = mcstate::random::binomial<real_type>(rng_state, I, p_IR);
    state_next[0] = S - n_SI;
    state_next[1] = I + n_SI - n_IR;
    state_next[2] = R + n_IR;
    state_next[3] = cases_cumul + n_SI;
    // state_next[4] = (time % shared.freq == 0) ? n_SI : (cases_inc + n_SI);
    state_next[4] = n_SI;
  }

  // Then, rather than a constructor we have some converters:
  static shared_state build_shared(cpp11::list pars) {
    const real_type I0 = with_default(10, pars["I0"]);
    const real_type N = with_default(1000, pars["N"]);

    const real_type beta = with_default(0.2, pars["beta"]);
    const real_type gamma = with_default(0.1, pars["gamma"]);

    const real_type exp_noise = with_default(1e6, pars["exp_noise"]);

    return shared_state{N, I0, beta, gamma, exp_noise};
  }

  // This one could be optional
  static internal_state build_internal(cpp11::list pars) {
    return sir::internal_state{};
  }

  static real_type compare_data(const real_type * state,
                                const data_type& data,
                                const shared_state& shared,
                                internal_state& internal,
                                rng_state_type& rng_state) {
    const auto incidence_observed = data.incidence;
    if (std::isnan(data.incidence)) {
      return 0;
    }
    const auto noise =
      mcstate::random::exponential(rng_state, shared.exp_noise);
    const auto incidence_modelled = state[4];
    const auto lambda = incidence_modelled + noise;
    return mcstate::density::poisson(incidence_observed, lambda, true);
  }
};
