#include <dust2/common.hpp>

// [[dust2::class(special)]]
// [[dust2::time_type(mixed)]]
// [[dust2::parameter(sd, required = FALSE, constant = FALSE)]]
class special {
public:
  special() = delete;

  using real_type = double;

  struct shared_state {
    real_type sd;
  };

  struct internal_state {};

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"a", {}}, {"b", {}}};
  }

  static size_t size_special() {
    return 1;
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    state_next[0] = 0;
    state_next[1] = 0;
  }

  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
    const real_type b = state[1];
    state_deriv[0] = b;
  }

  static void update(real_type time,
                     real_type dt,
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     rng_state_type& rng_state,
                     real_type * state_next) {
    const real_type y = state[1];
    state_next[1] = monty::random::normal<real_type>(rng_state, y, shared.sd);
  }

  static shared_state build_shared(cpp11::list pars) {
    const real_type sd = dust2::r::read_real(pars, "sd", 1.0);

    return shared_state{sd};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
  }
};
