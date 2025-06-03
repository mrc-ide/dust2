#include <dust2/common.hpp>

// [[dust2::class(specialevent)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(r0, required = FALSE, constant = TRUE)]]
// [[dust2::parameter(r1, required = FALSE, constant = TRUE)]]
// [[dust2::parameter(time_activate, required = FALSE, constant = TRUE)]]
class specialevent {
public:
  specialevent() = delete;

  using real_type = double;

  struct shared_state {
    real_type r0;
    real_type r1;
    real_type time_activate;
    bool use_time_activate;
  };

  struct internal_state {};

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"y", {}}, {"active", {}}, {"end", {}}};
  }

  static size_t size_special() {
    return 2;
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    state_next[0] = 0; // y
    state_next[1] = 0; // active
    state_next[2] = 0; // end
  }

  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
    const bool is_active = state[1] != 0;
    state_deriv[0] = is_active ? shared.r1 : shared.r0;
  }

  static auto events(const shared_state& shared, internal_state& internal) {
    dust2::ode::events_type<real_type> events;

    // There are a few events to consider here:
    //
    // * turn on at a time
    // * turn off after a duration
    // * turn on as we pass a point (up or down)
    //
    // Start this with just one that starts at a time to activate and
    // we can expand options later.
    if (shared.use_time_activate) {
      auto test = [&](const double t, const real_type* y) {
        return t - shared.time_activate;
      };
      auto action = [&](const double t, const double sign, double* y) {
        y[1] = 1;
      };
      events.push_back({"time_activate", {}, test, action});
    }
    return events;
  }

  static shared_state build_shared(cpp11::list pars) {
    const auto r0 = dust2::r::read_real(pars, "r0", 1);
    const auto r1 = dust2::r::read_real(pars, "r1", -1);

    const auto time_activate = dust2::r::read_real(pars, "time_activate", NA_REAL);
    const bool use_time_activate = !std::isnan(time_activate);

    return shared_state{r0, r1, time_activate, use_time_activate};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
  }
};
