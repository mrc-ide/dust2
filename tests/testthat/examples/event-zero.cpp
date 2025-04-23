#include <dust2/common.hpp>

// [[dust2::class(eventzero)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(r, rank = 0)]]
// [[dust2::parameter(value_y, rank = 0)]]
// [[dust2::parameter(value_z, rank = 0)]]
// [[dust2::parameter(period_z, rank = 0)]]
class eventzero {
public:
  eventzero() = delete;

  using real_type = double;

  struct shared_state {
    real_type r;
    real_type value_y;
    real_type value_z;
    real_type period_z;
 };

  struct internal_state {};

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"y", {}}, {"z", {}}};
  }

  static auto zero_every(const shared_state& shared) {
    // Reset 'z' every period_z, but we'l have 'y' accumulate
    return dust2::zero_every_type<real_type>{{shared.period_z, {1}}};
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state) {
    state[0] = 0;
    state[1] = 1;
  }

  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
    state_deriv[0] = time / shared.r;
    state_deriv[1] = time / shared.r;
  }

  static shared_state build_shared(cpp11::list pars) {
    const real_type r = dust2::r::read_real(pars, "r", 1.0);
    const real_type value_y = dust2::r::read_real(pars, "value_y", 9999.0);
    const real_type value_z = dust2::r::read_real(pars, "value_z", 9999.0);
    const real_type period_z = dust2::r::read_real(pars, "period_z", 1.0);
    return shared_state{r, value_y, value_z, period_z};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
  }

  static auto events(const shared_state& shared, internal_state& internal) {
    const auto test_y = [&](real_type t, const real_type* y) { return y[0] - shared.value_y; };
    const auto test_z = [&](real_type t, const real_type* y) { return y[0] - shared.value_z; };
    auto no_action = [](const double t, const double sign, double* y) {
    };

    dust2::ode::events_type<real_type> events;
    events.push_back(dust2::ode::event<real_type>("y", {0}, test_y, no_action, dust2::ode::root_type::increase));
    events.push_back(dust2::ode::event<real_type>("z", {1}, test_z, no_action, dust2::ode::root_type::increase));
    return events;
  }
};
