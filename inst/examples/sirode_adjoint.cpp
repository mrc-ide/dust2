#include <dust2/common.hpp>

// Continuous-time SIR with hand-written adjoint methods for testing
// the continuous adjoint backward integrator.
//
// States: S, I, R, cases_cumul, cases_inc (5 ODE variables)
// Parameters differentiated: beta, gamma, I0 (3 adjoint params)
// Adjoint vector: [adj_S, adj_I, adj_R, adj_cases_cumul, adj_cases_inc,
//                  adj_beta, adj_gamma, adj_I0]  (8 elements)
//
// The RHS is:
//   dS/dt         = -beta * S * I / N
//   dI/dt         = beta * S * I / N - gamma * I
//   dR/dt         = gamma * I
//   dcases_cumul/dt = beta * S * I / N
//   dcases_inc/dt   = beta * S * I / N
//
// The adjoint RHS computes (∂f/∂y)^T * λ_state + (∂f/∂θ)^T * λ_state
// for backward integration.

// [[dust2::class(sirode_adjoint)]]
// [[dust2::time_type(continuous)]]
// [[dust2::has_compare()]]
// [[dust2::has_adjoint()]]
// [[dust2::parameter(I0, rank = 0, constant = FALSE, required = FALSE)]]
// [[dust2::parameter(N, rank = 0, constant = TRUE, required = FALSE)]]
// [[dust2::parameter(beta, rank = 0, constant = FALSE, required = FALSE)]]
// [[dust2::parameter(gamma, rank = 0, constant = FALSE, required = FALSE)]]
// [[dust2::parameter(exp_noise, rank = 0, constant = TRUE, required = FALSE)]]
class sirode_adjoint {
public:
  sirode_adjoint() = delete;

  using real_type = double;

  struct shared_state {
    real_type N;
    real_type I0;
    real_type beta;
    real_type gamma;
    real_type exp_noise;
    struct {
      struct {
        dust2::packing state;
        dust2::packing adjoint;
      } packing;
    } odin;
  };

  struct internal_state {};

  struct data_type {
    real_type incidence;
  };

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"S", {}}, {"I", {}}, {"R", {}},
                          {"cases_cumul", {}}, {"cases_inc", {}}};
  }

  static dust2::packing packing_gradient(const shared_state& shared) {
    return dust2::packing{{"beta", {}}, {"gamma", {}}, {"I0", {}}};
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

  // Adjoint of the RHS: computes the backward-time derivatives of the
  // adjoint vector.  For the augmented backward system, this is:
  //
  //   d(adj_state)/dτ = (∂f/∂y)^T * adj_state
  //   d(adj_param)/dτ = (∂f/∂θ)^T * adj_state
  //
  // where f = [dS/dt, dI/dt, dR/dt, dcases_cumul/dt, dcases_inc/dt].
  //
  // Jacobian ∂f/∂y (rows = equations, cols = states):
  //         S                I              R  cases_cumul  cases_inc
  //  dS:  -beta*I/N       -beta*S/N        0       0           0
  //  dI:   beta*I/N   beta*S/N - gamma     0       0           0
  //  dR:      0            gamma            0       0           0
  //  cc:   beta*I/N        beta*S/N         0       0           0
  //  ci:   beta*I/N        beta*S/N         0       0           0
  //
  // (∂f/∂y)^T * λ means column j of J^T dotted with λ, i.e., row j of
  // J dotted with λ... no: (J^T)_{j,k} = J_{k,j}, so (J^T λ)_j = Σ_k J_{k,j} λ_k.
  // That is, column j of J (all equations' partial wrt state j) dotted with λ.
  //
  // Jacobian ∂f/∂θ (rows = equations, cols = params):
  //         beta         gamma      I0
  //  dS:  -S*I/N          0         0
  //  dI:   S*I/N         -I         0
  //  dR:      0            I         0
  //  cc:   S*I/N           0         0
  //  ci:   S*I/N           0         0
  static void adjoint_rhs(real_type time,
                          const real_type * state,
                          const real_type * adjoint,
                          const shared_state& shared,
                          internal_state& internal,
                          real_type * adjoint_deriv) {
    const auto S = state[0];
    const auto I = state[1];
    const auto adj_S = adjoint[0];
    const auto adj_I = adjoint[1];
    // adj_R = adjoint[2] — not needed, R column of J is all zeros
    // except gamma in row dR, handled below
    const auto adj_R = adjoint[2];
    const auto adj_cc = adjoint[3];
    const auto adj_ci = adjoint[4];

    const auto bSN = shared.beta * S / shared.N;
    const auto bIN = shared.beta * I / shared.N;
    const auto SIN = S * I / shared.N;

    // (∂f/∂y)^T * λ_state:
    // adj_deriv[0] = d/dS column: -bIN*adj_S + bIN*adj_I + bIN*adj_cc + bIN*adj_ci
    adjoint_deriv[0] = bIN * (-adj_S + adj_I + adj_cc + adj_ci);
    // adj_deriv[1] = d/dI column: -bSN*adj_S + (bSN-gamma)*adj_I + gamma*adj_R + bSN*adj_cc + bSN*adj_ci
    adjoint_deriv[1] = bSN * (-adj_S + adj_I + adj_cc + adj_ci) -
                       shared.gamma * adj_I + shared.gamma * adj_R;
    // adj_deriv[2] = d/dR column: all zeros
    adjoint_deriv[2] = 0;
    // adj_deriv[3] = d/d(cases_cumul) column: all zeros
    adjoint_deriv[3] = 0;
    // adj_deriv[4] = d/d(cases_inc) column: all zeros
    adjoint_deriv[4] = 0;

    // (∂f/∂θ)^T * λ_state:
    // adj_deriv[5] = d/d(beta): -SIN*adj_S + SIN*adj_I + SIN*adj_cc + SIN*adj_ci
    adjoint_deriv[5] = SIN * (-adj_S + adj_I + adj_cc + adj_ci);
    // adj_deriv[6] = d/d(gamma): -I*adj_I + I*adj_R
    adjoint_deriv[6] = I * (-adj_I + adj_R);
    // adj_deriv[7] = d/d(I0): 0 (I0 doesn't appear in RHS)
    adjoint_deriv[7] = 0;
  }

  static void adjoint_compare_data(real_type time,
                                   const real_type * state,
                                   const real_type * adjoint,
                                   const data_type& data,
                                   const shared_state& shared,
                                   internal_state& internal,
                                   real_type * adjoint_next) {
    // Copy adjoint through
    const size_t n = 8;
    for (size_t k = 0; k < n; ++k) {
      adjoint_next[k] = adjoint[k];
    }
    if (std::isnan(data.incidence)) {
      return;
    }
    const auto cases_inc = state[4];
    const auto noise = 1.0 / shared.exp_noise;
    const auto lambda = cases_inc + noise;
    // d/d(cases_inc) of Poisson(incidence | lambda):
    //   d/d(cases_inc) [incidence * log(lambda) - lambda]
    //   = incidence / lambda - 1
    const auto adj_lambda = data.incidence / lambda - 1;
    adjoint_next[4] += adj_lambda;
  }

  static void adjoint_initial(real_type time,
                              const real_type * state,
                              const real_type * adjoint,
                              const shared_state& shared,
                              internal_state& internal,
                              real_type * adjoint_next) {
    for (size_t k = 0; k < 8; ++k) {
      adjoint_next[k] = adjoint[k];
    }
    // initial(I) = I0, so d(initial_I)/d(I0) = 1
    // initial(S) = N - I0, so d(initial_S)/d(I0) = -1
    // adj_I0 += adj_I * 1 + adj_S * (-1)
    adjoint_next[7] += adjoint[1] - adjoint[0];
  }

  static shared_state build_shared(cpp11::list pars) {
    const real_type I0 = dust2::r::read_real(pars, "I0", 10);
    const real_type N = dust2::r::read_real(pars, "N", 1000);
    const real_type beta = dust2::r::read_real(pars, "beta", 0.2);
    const real_type gamma = dust2::r::read_real(pars, "gamma", 0.1);
    const real_type exp_noise = dust2::r::read_real(pars, "exp_noise", 1e6);
    auto shared = shared_state{N, I0, beta, gamma, exp_noise, {}};
    shared.odin.packing.state = packing_state(shared);
    shared.odin.packing.adjoint = dust2::packing{
      {"S", {}}, {"I", {}}, {"R", {}},
      {"cases_cumul", {}}, {"cases_inc", {}},
      {"beta", {}}, {"gamma", {}}, {"I0", {}}};
    return shared;
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.I0 = dust2::r::read_real(pars, "I0", shared.I0);
    shared.beta = dust2::r::read_real(pars, "beta", shared.beta);
    shared.gamma = dust2::r::read_real(pars, "gamma", shared.gamma);
  }

  static auto zero_every(const shared_state& shared) {
    return dust2::zero_every_type<real_type>{{1, {4}}};
  }

  static data_type build_data(cpp11::list r_data, const shared_state& shared) {
    auto data = static_cast<cpp11::list>(r_data);
    auto incidence = dust2::r::read_real(data, "incidence", NA_REAL);
    return data_type{incidence};
  }

  static real_type compare_data(const real_type time,
                                const real_type * state,
                                const data_type& data,
                                const shared_state& shared,
                                internal_state& internal,
                                rng_state_type& rng_state) {
    if (std::isnan(data.incidence)) {
      return 0;
    }
    const auto cases_inc = state[4];
    const auto noise = 1.0 / shared.exp_noise;
    const auto lambda = cases_inc + noise;
    return monty::density::poisson(data.incidence, lambda, true);
  }
};
