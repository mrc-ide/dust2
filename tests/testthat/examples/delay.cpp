#include <dust2/common.hpp>

// [[dust2::class(sirdelay)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(I0, rank = 0, constant = FALSE, required = FALSE)]]
// [[dust2::parameter(N, rank = 0, constant = TRUE, required = FALSE)]]
// [[dust2::parameter(beta, rank = 0, constant = FALSE, required = FALSE)]]
// [[dust2::parameter(gamma, rank = 0, constant = FALSE, required = FALSE)]]
class sirdelay {
public:
  sirdelay() = delete;

  using real_type = double;

  struct shared_state {
    real_type N;
    real_type I0;
    real_type beta;
    real_type gamma;
  };

  struct internal_state {};

  struct data_type {
    real_type incidence;
  };

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"S", {}}, {"I", {}}, {"R", {}}, {"cases_cumul", {}}, {"cases_inc", {}}};
  }

  static size_t size_output() {
    return 1;
  }

  static auto delays(const shared_state& shared) {
    return dust2::ode::delays<real_type>(false, true, {{1, {3}}});
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
                  const dust2::ode::delay_result_type<real_type>& delays,
                  real_type * state_deriv) {
    const auto S = state[0];
    const auto I = state[1];
    const auto rate_SI = shared.beta * S * I / shared.N;
    const auto rate_IR = shared.gamma * I;
    state_deriv[0] = -rate_SI;
    state_deriv[1] = rate_SI - rate_IR;
    state_deriv[2] = rate_IR;
    state_deriv[3] = rate_SI;
  }

  static void output(real_type time,
                     real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     const dust2::ode::delay_result_type<real_type>& delays) {
    state[4] = state[3] - delays[0][0];
  }

  static shared_state build_shared(cpp11::list pars) {
    const real_type I0 = dust2::r::read_real(pars, "I0", 10);
    const real_type N = dust2::r::read_real(pars, "N", 1000);
    const real_type beta = dust2::r::read_real(pars, "beta", 0.2);
    const real_type gamma = dust2::r::read_real(pars, "gamma", 0.1);
    return shared_state{N, I0, beta, gamma};
  }

  // This is the bit that we'll use to do fast parameter updating, and
  // we'll guarantee somewhere that the size does not change.
  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.I0 = dust2::r::read_real(pars, "I0", shared.I0);
    shared.beta = dust2::r::read_real(pars, "beta", shared.beta);
    shared.gamma = dust2::r::read_real(pars, "gamma", shared.gamma);
  }
};
