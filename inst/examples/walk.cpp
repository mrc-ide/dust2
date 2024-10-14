#include <dust2/common.hpp>

// [[dust2::class(walk)]]
// [[dust2::time_type(discrete)]]
// [[dust2::parameter(sd)]]
// [[dust2::parameter(len)]]
// [[dust2::parameter(random_initial)]]
class walk {
public:
  // No constructor - turning this off is optional
  walk() = delete;

  // We need to be able to swap between different precisions when using
  using real_type = double;

  // These are parameters that do not vary within a block.  There's a
  // good question here of if we need a further layer, which would
  // help with running simulations across many parameters or on a GPU.
  //
  // This will be decided pretty soon.  The big issue is that we'll
  // have quite a bit of variation in exactly what is wanted to be
  // changed in different contexts.
  struct shared_state {
    size_t len;
    real_type sd;
    bool random_initial;
  };

  // This would be data that is specific to a particle.  There is no
  // guarantee that a particle will see the same copy of this in
  // sequential calls.  We might want to call this 'scratch' rather
  // than internal really.
  struct internal_state {};

  // This one always feels a bit weird, really.
  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"x", {shared.len}}};
  }

  // This is the bit that we'll use to do fast parameter updating, and
  // we'll guarantee somewhere that the size does not change.
  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.sd = dust2::r::read_real(pars, "sd", shared.sd);
  }

  // Compute initial state.  It's not completely clear that we would
  // want dt here, but I am including it for now at least.
  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    if (shared.random_initial) {
      for (size_t i = 0; i < shared.len; ++i) {
        state_next[i] =
          monty::random::normal<real_type>(rng_state, 0, shared.sd);
      }
    } else {
      std::fill(state_next, state_next + shared.len, 0);
    }
  }

  // The main update function, converting state to state_next
  static void update(real_type time,
                     real_type dt,
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     rng_state_type& rng_state,
                     real_type * state_next) {
    for (size_t i = 0; i < shared.len; ++i) {
      state_next[i] =
        monty::random::normal(rng_state, state[i], shared.sd * dt);
    }
  }

  // Then, rather than a constructor we have some converters:
  static shared_state build_shared(cpp11::list pars) {
    const auto len = dust2::r::read_size(pars, "len", 1);
    const auto sd = dust2::r::read_real(pars, "sd");
    const auto random_initial =
      dust2::r::read_bool(pars, "random_initial", false);
    return shared_state{len, sd, random_initial};
  }
};
