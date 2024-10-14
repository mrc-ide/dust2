#include <dust2/common.hpp>

// [[dust2::class(walk)]]
// [[dust2::time_type(discrete)]]
// [[dust2::parameter(sd)]]
// [[dust2::parameter(len)]]
class walk {
public:
  walk() = delete;

  using real_type = double;
  using rng_state_type = monty::random::generator<real_type>;

  struct shared_state {
    size_t len;
    real_type sd;
  };

  struct internal_state {
    std::vector<real_type> scratch;
  };

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"x", {}}};
  }

  static shared_state build_shared(cpp11::list pars) {
    const auto len = dust2::r::read_size(pars, "len", 10);
    const auto sd = dust2::r::read_real(pars, "sd");
    return shared_state{len, sd};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.sd = dust2::r::read_real(pars, "sd", shared.sd);
  }

  static internal_state build_internal(const shared_state& shared) {
    std::vector<real_type> scratch(shared.len);
    return internal_state{scratch};
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    std::fill(state_next, state_next + shared.len, 0);
  }

  static void update(real_type time,
                     real_type dt,
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     rng_state_type& rng_state,
                     real_type * state_next) {
    const auto x = state[0];
    for (size_t i = 0; i < shared.len; ++i) {
      internal.scratch[i] =
        monty::random::normal(rng_state, x, shared.sd * dt);
    }
    state_next[0] = std::accumulate(internal.scratch.begin(),
                                    internal.scratch.end(),
                                    static_cast<real_type>(0.0)) / shared.len;
  }
};
