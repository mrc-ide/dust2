#include <dust2/common.hpp>

// [[dust2::class(change)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(r, rank = 0)]]
// [[dust2::parameter(n, rank = 0)]]
// [[dust2::parameter(t_change, rank = 1)]]
// [[dust2::parameter(delta, rank = 1)]]
class change {
public:
  change() = delete;

  using real_type = double;

  struct shared_state {
    real_type r;
    std::vector<real_type> t_change;
    std::vector<real_type> delta;
  };

  struct internal_state {};

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"y", {}}};
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state) {
    state[0] = 0;
  }

  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
    state_deriv[0] = shared.r;
  }

  static shared_state build_shared(cpp11::list pars) {
    const real_type r1 = dust2::r::read_real(pars, "r");
    const size_t n = dust2::r::read_int(pars, "n");
    std::vector<real_type> t_change(n);
    std::vector<real_type> delta(n);
    dust2::r::read_real_vector(pars, n, t_change.data(), "t_change", true);
    dust2::r::read_real_vector(pars, n, delta.data(), "delta", false);
    return shared_state{r1, t_change, delta};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
  }

  static auto events(const shared_state& shared, internal_state& internal) {
    const auto n = shared.t_change.size();
    dust2::ode::events_type<real_type> events;
    events.reserve(n);
    std::string prefix = "event-";
    for (size_t i = 0; i < n; ++i) {
      const auto name = prefix + std::to_string(i + 1);
      auto test = [&, i](const double t, const real_type* y) {
        return t - shared.t_change[i];
      };
      auto action = [&, i](const double t, const double sign, double* y) {
        y[0] += shared.delta[i];
      };
      events.push_back(dust2::ode::event<real_type>(name, {}, test, action));
    }
    return events;
  }
};
