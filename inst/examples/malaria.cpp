#include <dust2/common.hpp>

// [[dust2::class(malaria)]]
// [[dust2::time_type(mixed)]]
// [[dust2::has_compare()]]
// [[dust2::parameter(a, required = FALSE, constant = FALSE)]]
// [[dust2::parameter(n_rates, required = FALSE, constant = TRUE, type = "int")]]
// [[dust2::parameter(r, required = FALSE, constant = FALSE)]]
// [[dust2::parameter(tau, required = FALSE, constant = FALSE)]]
// [[dust2::parameter(initial_Ih, required = FALSE, constant = TRUE)]]
// [[dust2::parameter(initial_Iv, required = FALSE, constant = TRUE)]]
// [[dust2::parameter(initial_Sv, required = FALSE, constant = TRUE)]]
class malaria {
public:
  malaria() = delete;

  using real_type = double;

  struct shared_state {
    real_type a; // Biting rate (bites per human per mosquito)

    real_type bh;   // Pr(transmission vector to human)
    real_type bv;   // Pr(transmission human to vector)
    size_t n_rates; // number of exposure compartments
    real_type mu;   // -log(Pr(vector survival))
    real_type r;    // Rate of recovery
    real_type tau;  // Length in mosquito latency period

    real_type initial_Ih; // Initial infected humans
    real_type initial_Iv; // Initial infected vectors
    real_type initial_Sv; // Initial susceptible vector

    real_type beta_volatility; // Volatility of random walk
  };

  struct internal_state {};

  struct data_type {
    real_type tested;
    real_type positive;
  };

  using rng_state_type = monty::random::generator<real_type>;

  // State is modelled as [Sh, Ih, Sv, Iv, beta, Ev...]
  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"Sh", {}}, {"Ih", {}}, {"Sv", {}}, {"Iv", {}}, {"beta", {}}, {"Ev", {shared.n_rates}}};
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    state_next[0] = 1 - shared.initial_Ih; // Sh
    state_next[1] = shared.initial_Ih;     // Ih
    state_next[2] = shared.initial_Sv;     // Sv
    state_next[3] = shared.initial_Iv;     // Iv
    state_next[4] = shared.mu;             // beta
    for (size_t i = 5; i < 5 + shared.n_rates; ++i) {
      state_next[i] = 0;
    }
  }

  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
    const real_type Sh = state[0];
    const real_type Ih = state[1];
    const real_type Sv = state[2];
    const real_type Iv = state[3];
    const real_type beta = state[4];
    const real_type * Ev = state + 5;
    const real_type foi_h = shared.a * shared.bh * Iv;
    const real_type foi_v = shared.a * shared.bv * Ih;
    state_deriv[5] = foi_v * Sv - (shared.n_rates / shared.tau) * Ev[0] - shared.mu * Ev[0];
    for (size_t i = 1; i < shared.n_rates; ++i) {
      state_deriv[5 + i] = (shared.n_rates / shared.tau) * Ev[i - 1] - (shared.n_rates / shared.tau) * Ev[i] - shared.mu * Ev[i];
    }
    state_deriv[0] = - foi_h * Sh + shared.r * Ih;
    state_deriv[1] = foi_h * Sh - shared.r * Ih;
    state_deriv[2] = beta * shared.initial_Sv - foi_v * Sv - shared.mu * Sv;
    state_deriv[3] = (shared.n_rates / shared.tau) * Ev[shared.n_rates - 1] - shared.mu * Iv;
  }

  static void update(real_type time,
                     real_type dt,
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     rng_state_type& rng_state,
                     real_type * state_next) {
    const real_type beta = state[4];
    state_next[4] = beta * monty::math::exp(monty::random::normal<real_type>(rng_state, 0, shared.beta_volatility));
  }

  static shared_state build_shared(cpp11::list pars) {
    const real_type a = dust2::r::read_real(pars, "a", 1.0 / 3.0);
    const size_t n_rates = dust2::r::read_size(pars, "n_rates", 15);
    const real_type r = dust2::r::read_real(pars, "r", 0.01);
    const real_type tau = dust2::r::read_real(pars, "tau", 12);

    const real_type initial_Ih = dust2::r::read_real(pars, "initial_Ih", 0.8);
    const real_type initial_Iv = dust2::r::read_real(pars, "initial_Iv", 18);
    const real_type initial_Sv = dust2::r::read_real(pars, "initial_Sv", 100);

    const real_type beta_volatility = 0.5;
    const real_type bh = 0.05;
    const real_type bv = 0.05;
    const real_type p = 0.9; // Daily probability of vector survival
    const real_type mu = - monty::math::log(p);

    return shared_state{a, bh, bv, n_rates, mu, r, tau, initial_Ih, initial_Iv, initial_Sv, beta_volatility};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.a = dust2::r::read_real(pars, "a", shared.a);
    shared.r = dust2::r::read_real(pars, "r", shared.r);
    shared.tau = dust2::r::read_real(pars, "tau", shared.tau);
  }

  static data_type build_data(cpp11::list r_data, const shared_state& shared) {
    auto data = static_cast<cpp11::list>(r_data);
    auto tested = dust2::r::read_real(data, "tested", NA_REAL);
    auto positive = dust2::r::read_real(data, "positive", NA_REAL);
    return data_type{tested, positive};
  }

  static real_type compare_data(const real_type time,
                                const real_type * state,
                                const data_type& data,
                                const shared_state& shared,
                                internal_state& internal,
                                rng_state_type& rng_state) {

    if (std::isnan(data.positive) || std::isnan(data.tested)) {
      return 0;
    }
    const real_type Ih = state[1];
    return monty::density::binomial(data.positive, data.tested, Ih, true);
  }
};
