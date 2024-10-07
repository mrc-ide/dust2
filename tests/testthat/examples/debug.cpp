#include <dust2/common.hpp>
#include <dust2/r/debug.hpp>

// [[dust2::class(sirdebug)]]
// [[dust2::time_type(discrete)]]
// [[dust2::parameter(I0)]]
// [[dust2::parameter(N)]]
// [[dust2::parameter(beta)]]
// [[dust2::parameter(gamma)]]
// [[dust2::parameter(exp_noise)]]
class sirdebug {
public:
  sirdebug() = delete;

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

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"S", {}}, {"I", {}}, {"R", {}}, {"cases_cumul", {}}, {"cases_inc", {}}};
  }

  static dust2::packing packing_gradient(const shared_state& shared) {
    return dust2::packing{};
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
    const auto cases_inc = state[4];
    const auto p_SI = 1 - monty::math::exp(-shared.beta * I / shared.N * dt);
    const auto p_IR = 1 - monty::math::exp(-shared.gamma * dt);
    const auto n_SI = monty::random::binomial<real_type>(rng_state, S, p_SI);
    const auto n_IR = monty::random::binomial<real_type>(rng_state, I, p_IR);
    state_next[0] = S - n_SI;
    state_next[1] = I + n_SI - n_IR;
    state_next[2] = R + n_IR;
    state_next[3] = cases_cumul + n_SI;
    state_next[4] = cases_inc + n_SI;

    if (time > 2 && time < 5) {
      auto env = dust2::r::debug::create_env();
      dust2::r::debug::save(time, "time", env);
      dust2::r::debug::save(S, "S", env);
      dust2::r::debug::save(I, "I", env);
      dust2::r::debug::save(R, "R", env);
      dust2::r::debug::save(n_SI, "n_SI", env);
      dust2::r::debug::save(n_IR, "n_IR", env);
      dust2::r::debug::browser(env, "update", time);
    }
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

  static auto zero_every(const shared_state& shared) {
    return dust2::zero_every_type<real_type>{{1, {4}}}; // zero[1] = {4};
  }
};
