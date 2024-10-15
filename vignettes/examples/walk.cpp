#include <dust2/common.hpp>

// [[dust2::class(walk)]]
// [[dust2::time_type(discrete)]]
// [[dust2::parameter(sd)]]
class walk {
public:
  walk() = delete;

  using real_type = double;

  struct shared_state {
    real_type sd;
  };

  struct internal_state {};

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"x", {}}};
  }

  static shared_state build_shared(cpp11::list pars) {
    const auto sd = dust2::r::read_real(pars, "sd");
    return shared_state{sd};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.sd = dust2::r::read_real(pars, "sd", shared.sd);
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    state_next[0] = 0;
  }

  static void update(real_type time,
                     real_type dt,
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     rng_state_type& rng_state,
                     real_type * state_next) {
    state_next[0] = monty::random::normal(rng_state, state[0], shared.sd * dt);
  }
};
