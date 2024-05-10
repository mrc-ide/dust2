#include <dust2/common.hpp>

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
  using internal_state = dust2::no_internal_state;

  // This model cannot be used with a particle filter; it has no
  // internal data type.
  using data_type = dust2::no_data;

  // This one always feels a bit weird, really.
  using rng_state_type = mcstate::random::generator<real_type>;

  static size_t size(const shared_state& shared) {
    return shared.len;
  }

  // This is the bit that we'll use to do fast parameter updating, and
  // we'll guarantee somewhere that the size does not change.
  static void update_shared(cpp11::list pars, shared_state& shared) {
    const cpp11::sexp r_sd = pars["sd"];
    if (r_sd != R_NilValue) {
      shared.sd = cpp11::as_cpp<walk::real_type>(pars["sd"]);
    }
  }

  // This is a reasonable default implementation in the no-internal
  // case
  static void update_internal(const shared_state& shared,
                              internal_state& internal) {
  }

  // Compute initial state.  It's not completely clear that we would
  // want dt here, but I am including it for now at least.
  static void initial(real_type time,
                      real_type dt,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    if (shared.random_initial) {
      for (size_t i = 0; i < shared.len; ++i) {
        state_next[i] =
          mcstate::random::normal<real_type>(rng_state, 0, shared.sd);
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
        mcstate::random::normal(rng_state, state[i], shared.sd * dt);
    }
  }

  // Then, rather than a constructor we have some converters:
  static shared_state build_shared(cpp11::list pars) {
    size_t len = 1;
    const cpp11::sexp r_len = pars["len"];
    if (r_len != R_NilValue) {
      len = cpp11::as_cpp<int>(r_len);
    }
    const walk::real_type sd = cpp11::as_cpp<walk::real_type>(pars["sd"]);
    const bool random_initial = pars["random_initial"] == R_NilValue ? false :
      cpp11::as_cpp<bool>(pars["random_initial"]);
    return shared_state{len, sd, random_initial};
  }

  // This one could be optional
  static internal_state build_internal(cpp11::list pars) {
    return walk::internal_state{};
  }
};
