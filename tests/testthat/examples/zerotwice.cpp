#include <dust2/common.hpp>

// [[dust2::class(zerotwice)]]
// [[dust2::time_type(discrete)]]
class zerotwice {
public:
  zerotwice() = delete;

  using real_type = double;

  struct shared_state {
    real_type a;
    real_type b;
  };
  using internal_state = dust2::no_internal_state;
  using data_type = dust2::no_data;
  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"x", {}}, {"y", {}}};
  }

  static dust2::packing packing_gradient(const shared_state& shared) {
    return dust2::packing{};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    state_next[0] = 0;
    state_next[1] = 0;
  }

  static void update(real_type time,
                     real_type dt,
                     const real_type * state,
                     const shared_state& shared,
                     internal_state& internal,
                     rng_state_type& rng_state,
                     real_type * state_next) {
    state_next[0] = state[0] + 1;
    state_next[1] = state[1] + 1;
  }

  static shared_state build_shared(cpp11::list pars) {
    const auto a = dust2::r::read_real(pars, "a");
    const auto b = dust2::r::read_real(pars, "b");
    return shared_state{a, b};
  }

  static internal_state build_internal(const shared_state& shared) {
    return internal_state{};
  }

  static auto zero_every(const shared_state& shared) {
    return dust2::zero_every_type<real_type>{{shared.a, {0}}, {shared.b, {1}}};
  }
};
