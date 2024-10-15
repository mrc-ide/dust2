#include <dust2/common.hpp>

// [[dust2::class(volatility)]]
// [[dust2::time_type(discrete)]]
// [[dust2::has_compare()]]
// [[dust2::parameter(alpha)]]
// [[dust2::parameter(sigma)]]
// [[dust2::parameter(gamma)]]
// [[dust2::parameter(tau)]]
class volatility {
public:
  volatility() = delete;

  using real_type = double;

  struct shared_state {
    real_type alpha;
    real_type sigma;
    real_type gamma;
    real_type tau;
  };

  struct internal_state {};

  struct data_type {
    real_type observed;
  };

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"x", {}}};
  }

  static shared_state build_shared(cpp11::list pars) {
    const auto alpha = dust2::r::read_real(pars, "alpha", 0.91);
    const auto sigma = dust2::r::read_real(pars, "sigma", 1);
    const auto gamma = dust2::r::read_real(pars, "gamma", 1);
    const auto tau = dust2::r::read_real(pars, "tau", 1);
    return shared_state{alpha, sigma, gamma, tau};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.alpha = dust2::r::read_real(pars, "alpha", shared.alpha);
    shared.sigma = dust2::r::read_real(pars, "sigma", shared.sigma);
    shared.gamma = dust2::r::read_real(pars, "gamma", shared.gamma);
    shared.tau = dust2::r::read_real(pars, "tau", shared.tau);
  }

  static data_type build_data(cpp11::list r_data, const shared_state& shared) {
    auto data = static_cast<cpp11::list>(r_data);
    const auto observed = dust2::r::read_real(data, "observed", NA_REAL);
    return data_type{observed};
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    state_next[0] = monty::random::normal<real_type>(rng_state, 0, 1);
  }

  static void update(real_type time,
                     real_type dt,
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     rng_state_type& rng_state,
                     real_type * state_next) {
    const auto x = state[0];
    state_next[0] = shared.alpha * x +
      shared.sigma * monty::random::normal<real_type>(rng_state, 0, 1);
  }

  static real_type compare_data(const real_type time,
                                const real_type * state,
                                const data_type& data,
                                const shared_state& shared,
                                internal_state& internal,
                                rng_state_type& rng_state) {
    const auto x = state[0];
    return monty::density::normal(data.observed, shared.gamma * x, shared.tau,
                                  true);
  }
};
