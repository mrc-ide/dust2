#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <lostturnip.hpp>

#include <dust2/tools.hpp>
#include <dust2/zero.hpp>
#include <dust2/continuous/control.hpp>
#include <dust2/continuous/events.hpp>
#include <dust2/continuous/history.hpp>

namespace dust2 {
namespace ode {

template <typename T>
T square(T x) {
  return x * x;
}

template <typename T>
T clamp(T x, T min, T max) {
  return std::max(std::min(x, max), min);
}

// This is the internal state separate from 'y' that defines the
// system.  This includes the derivatives.
//
// Pulling this out into a little struct means that we can copy state
// from one particle to another fairly efficiently
template <typename real_type>
struct internals {
  history_step<real_type> last;
  history<real_type> history_values;
  event_history<real_type> events;

  std::vector<real_type> dydt;
  std::vector<real_type> step_times;
  real_type step_size;
  real_type error;
  size_t n_steps;
  size_t n_steps_accepted;
  size_t n_steps_rejected;
  bool save_history_;

  internals(size_t n_variables, size_t n_special, bool save_history) :
    last(n_variables, n_special),
    history_values(n_variables + n_special),
    dydt(n_variables),
    save_history_(save_history) {
    reset(last.c1.data());
  }

  void save_history() {
    if (save_history_) {
      history_values.add(last);
    }
  }

  void reset(const real_type * y) {
    step_size = 0;
    error = 0;
    n_steps = 0;
    n_steps_accepted = 0;
    n_steps_rejected = 0;
    step_times.clear();
    history_values.reset(y);
  }
};

template <typename real_type>
class solver {
public:
  ode::control<real_type> control;

  solver(size_t n_variables, size_t n_special, ode::control<real_type> control) :
    control(control),
    n_variables_(n_variables),
    n_special_(n_special),
    y_next_(n_variables_ + n_special_),
    y_stiff_(n_variables_ + n_special_),
    k2_(n_variables_),
    k3_(n_variables_),
    k4_(n_variables_),
    k5_(n_variables_),
    k6_(n_variables_),
    facc1_(1 / control.factor_min),
    facc2_(1 / control.factor_max) {
  }

  template <typename Rhs>
  real_type try_step(const real_type t, const real_type h, const real_type* y,
                     const real_type* dydt, real_type* c5, Rhs rhs) {
    const real_type* k1 = dydt;
    const auto t_next = t + h;
    for (size_t i = 0; i < n_variables_; ++i) { // 22
      y_next_[i] = y[i] + h * A21 * k1[i];
    }
    rhs(t + C2 * h, y_next_.data(), k2_.data());
    for (size_t i = 0; i < n_variables_; ++i) { // 23
      y_next_[i] = y[i] + h * (A31 * k1[i] + A32 * k2_[i]);
    }
    rhs(t + C3 * h, y_next_.data(), k3_.data());
    for (size_t i = 0; i < n_variables_; ++i) { // 24
      y_next_[i] = y[i] + h * (A41 * k1[i] + A42 * k2_[i] + A43 * k3_[i]);
    }
    rhs(t + C4 * h, y_next_.data(), k4_.data());
    for (size_t i = 0; i < n_variables_; ++i) { // 25
      y_next_[i] = y[i] + h * (A51 * k1[i] + A52 * k2_[i] + A53 * k3_[i] +
                              A54 * k4_[i]);
    }
    rhs(t + C5 * h, y_next_.data(), k5_.data());
    for (size_t i = 0; i < n_variables_; ++i) { // 26
      y_stiff_[i] = y[i] + h * (A61 * k1[i]  + A62 * k2_[i] +
                               A63 * k3_[i] + A64 * k4_[i] +
                               A65 * k5_[i]);
    }
    rhs(t_next, y_stiff_.data(), k6_.data());
    for (size_t i = 0; i < n_variables_; ++i) { // 27
      y_next_[i] = y[i] + h * (A71 * k1[i]  + A73 * k3_[i] + A74 * k4_[i] +
                              A75 * k5_[i] + A76 * k6_[i]);
    }
    rhs(t_next, y_next_.data(), k2_.data());

    for (size_t i = 0; i < n_variables_; ++i) {
      c5[i] = h * (D1 * k1[i]  + D3 * k3_[i] + D4 * k4_[i] +
                   D5 * k5_[i] + D6 * k6_[i] + D7 * k2_[i]);
    }

    for (size_t i = 0; i < n_variables_; ++i) {
      k4_[i] = h * (E1 * k1[i]  + E3 * k3_[i] + E4 * k4_[i] +
                    E5 * k5_[i] + E6 * k6_[i] + E7 * k2_[i]);
    }

    // Compute error:
    const auto atol = control.atol;
    const auto rtol = control.rtol;
    real_type err = 0.0;
    for (size_t i = 0; i < n_variables_; ++i) {
      auto sk =
        k4_[i] / (atol + rtol * std::max(std::abs(y[i]), std::abs(y_next_[i])));
      err += square(sk);
    }
    return std::sqrt(err / n_variables_);
  }

  // Take a single step
  template <typename Rhs>
  real_type step(real_type t, real_type t_end, real_type* y,
                 const events_type<real_type>& events,
                 ode::internals<real_type>& internals, Rhs rhs) {
    auto success = false;
    auto reject = false;
    auto truncated = false;
    auto event = false;
    auto h = internals.step_size;

    while (!success) {
      if (internals.n_steps > control.max_steps) {
        // throw a nicer error for all of these, with the current
        // time etc.
        throw std::runtime_error("too many steps");
      }
      if (h < control.step_size_min) {
        throw std::runtime_error("step too small");
      }
      if (h <= std::abs(t) * std::numeric_limits<real_type>::epsilon()) {
        throw std::runtime_error("step size vanished");
      }
      truncated = t + h > t_end;
      if (truncated) {
        h = t_end - t;
      }

      const auto err = try_step(t, h, y, internals.dydt.data(),
                                internals.last.c5.data(), rhs);
      internals.n_steps++;
      const auto fac11 = std::pow(err, control.constant);

      if (err <= 1) {
        success = true;
        update_interpolation(t, h, y, internals);
        if (!events.empty()) {
          // Uses C++17's structured binding, which is a pecular syntax...
          // https://en.cppreference.com/w/cpp/language/structured_binding.html
          const auto [had_event, t_next] = apply_events(t, h, y, events, internals);
          if (had_event) {
            event = true;
            truncated = false;
            h = t_next - t;
            if (n_special_ > 0) {
              const auto src = y_next_.begin() + n_variables_;
              std::copy_n(src, n_special_, y + n_variables_);
              std::copy_n(src, n_special_, y_stiff_.begin() + n_variables_);
            }
            rhs(t_next, y_next_.data(), k2_.data());
          }
        }
        accept(t, h, y, internals);
        internals.n_steps_accepted++;
        if (control.debug_record_step_times) {
          internals.step_times.push_back(truncated ? t_end : t + h);
        }
        internals.save_history();
        if (!truncated && !event) {
          const auto fac_old =
            std::max(internals.error, static_cast<real_type>(1e-4));
          auto fac = fac11 / std::pow(fac_old, control.beta);
          fac = clamp(fac / control.factor_safe, facc2_, facc1_);
          const auto h_new = h / fac;
          const auto h_max = reject ? h : control.step_size_max;
          internals.step_size = std::min(h_new, h_max);
          internals.error = err;
        }
      } else {
        reject = true;
        if (internals.n_steps_accepted >= 1) {
          internals.n_steps_rejected++;
        }
        h /= std::min(facc1_, fac11 / control.factor_safe);
      }
    }

    return truncated ? t_end : t + h;
  }

  template <typename Rhs>
  void run(real_type t, real_type t_end, real_type* y,
           zero_every_type<real_type>& zero_every,
           const events_type<real_type>& events,
           ode::internals<real_type>& internals, Rhs rhs) {
    if (n_special_ > 0) {
      std::copy_n(y + n_variables_, n_special_, y_next_.begin() + n_variables_);
      std::copy_n(y + n_variables_, n_special_, y_stiff_.begin() + n_variables_);
    }
    if (control.critical_times.empty()) {
      while (t < t_end) {
        apply_zero_every(t, y, zero_every, internals);
        t = step(t, t_end, y, events, internals, rhs);
      }
    } else {
      // Slightly more complex loop which ensures we never integrate
      // over the times within our critical times.  The upper loop is
      // a special case of this but is kept simple.
      auto tc_end = control.critical_times.end();
      auto tc = std::upper_bound(control.critical_times.begin(), tc_end, t);
      auto t_end_i = (tc == tc_end || *tc >= t_end) ? t_end : *tc;
      while (t < t_end) {
        apply_zero_every(t, y, zero_every, internals);
        t = step(t, t_end_i, y, events, internals, rhs);
        if (t >= t_end_i && t < t_end) {
          ++tc;
          t_end_i = (tc == tc_end || *tc >= t_end) ? t_end : *tc;
        }
      }
    }
    apply_zero_every_final(t, y, zero_every, internals);
  }

  template <typename Rhs>
  void initialise(const real_type t, const real_type* y,
                  ode::internals<real_type>& internals, Rhs rhs) {
    internals.reset(y);
    if (control.debug_record_step_times) {
      internals.step_times.push_back(t);
    }
    auto f0 = internals.dydt.data();
    auto f1 = k3_.data();
    auto y1 = y_next_.data();
    if (n_special_ > 0) {
      std::copy_n(y + n_variables_, n_special_, y1 + n_variables_);
    }

    // Compute a first guess for explicit Euler as
    //   h = 0.01 * norm (y0) / norm (f0)
    // the increment for explicit euler is small compared to the solution
    rhs(t, y, f0);

    real_type norm_f = 0.0;
    real_type norm_y = 0.0;

    for (size_t i = 0; i < n_variables_; ++i) {
      const real_type sk = control.atol + control.rtol * std::abs(y[i]);
      norm_f += square(f0[i] / sk);
      norm_y += square(y[i] / sk);
    }
    // Magic numbers here, from Hairer
    real_type h = (norm_f <= 1e-10 || norm_y <= 1e-10) ?
      1e-6 : std::sqrt(norm_y / norm_f) * 0.01;
    h = std::min(h, control.step_size_max);

    // Perform an explicit Euler step
    for (size_t i = 0; i < n_variables_; ++i) {
      y1[i] = y[i] + h * f0[i];
    }
    rhs(t + h, y1, f1);

    // Estimate the second derivative of the solution:
    real_type der2 = 0.0;
    for (size_t i = 0; i < n_variables_; ++i) {
      const real_type sk = control.atol + control.rtol * std::abs(y[i]);
      der2 += square((f1[i] - f0[i]) / sk);
    }
    der2 = std::sqrt(der2) / h;

    // Step size is computed such that
    //   h^order * max(norm(f0), norm(der2)) = 0.01
    constexpr real_type order = 5;
    const real_type der12 = std::max(std::abs(der2), std::sqrt(norm_f));
    const real_type h1 = (der12 <= 1e-15) ?
      std::max(1e-6, std::abs(h) * 1e-3) :
      std::pow(0.01 / der12, 1.0 / order);
    h = std::min(100 * std::abs(h), h1);
    if (!std::isfinite(h)) {
      throw std::runtime_error("Initial step size was not finite");
    }
    h = std::max(h, control.step_size_min * 100);
    h = std::min(h, control.step_size_max);
    internals.step_size = h;
  }

private:
  void update_interpolation(real_type t, real_type h, real_type* y, ode::internals<real_type>& internals) {
    internals.last.t0 = t;
    internals.last.t1 = t + h;
    internals.last.h = h;
    for (size_t i = 0; i < n_variables_; ++i) {
      const auto ydiff = y_next_[i] - y[i];
      const auto bspl = h * internals.dydt[i] - ydiff;
      internals.last.c1[i] = y[i];
      internals.last.c2[i] = ydiff;
      internals.last.c3[i] = bspl;
      internals.last.c4[i] = -h * k2_[i] + ydiff - bspl;
    }
    // If there are any specials, we need to copy these over too.
    if (n_special_ > 0) {
      std::copy_n(y + n_variables_, n_special_,
                  internals.last.c1.begin() + n_variables_);
    }
  }

  void accept(real_type t, real_type h, real_type* y, ode::internals<real_type>& internals) {
    std::copy_n(k2_.begin(), n_variables_, internals.dydt.begin());
    std::copy_n(y_next_.begin(), n_variables_, y);
  }

  std::tuple<bool, real_type> apply_events(real_type t0, real_type h, const real_type* y,
                                           const events_type<real_type>& events,
                                           ode::internals<real_type>& internals) {
    real_type t1 = t0 + h;

    // It might be worth saving this storage space in the solver, but
    // I doubt it matters in pratice.  We need to save all the found
    // events and their signs, even though probably only one will be
    // found, because of the possibility that multiple events could
    // happen at the same time (e.g., two time-scheduled events or two
    // functions that happen to hit their roots at the same time).
    std::vector<bool> found(events.size());
    std::vector<real_type> sign(events.size());
    bool found_any = false;

    for (size_t idx_event = 0; idx_event < events.size(); ++idx_event) {
      const auto& e = events[idx_event];
      // Use y_stiff as temporary space here, it's only used
      // transiently and within the step
      real_type * y_t = y_stiff_.data();
      auto fn = [&](auto t) {
        internals.last.interpolate(t, e.index, y_t);
        return e.test(t, y_t);
      };
      const auto f_t0 = fn(t0);
      const auto f_t1 = fn(t1);
      if (is_root(f_t0, f_t1, e.root)) {
        // These probably should move into the ode control, but there
        // should really be any great need to change them, and the
        // interpolation is expected to be quite fast and accurate.
        constexpr real_type eps = 1e-6;
        constexpr size_t steps = 100;
        t1 = lostturnip::find_result<real_type>(fn, t0, t1, eps, steps).x;
        // Currently untested - in the case where we have two roots
        // that would have been crossed in this time window, the one
        // we are currently considering happens first so pre-empts the
        // previously found events.
        if (found_any) {
          std::fill(found.begin(), found.end(), false);
        }
      } else if (!(f_t1 == 0 && f_t0 != 0)) {
        // Consider the case where jump to a root *exactly* at t1;
        // this happens in coincident roots and with roots that are
        // based in time, and which we arrange for the solver to stop
        // at (e.g., while using simulate()).
        //
        // This test is the inverse of this though, because in the
        // case where we *don't* get an exact root we should skip the
        // bookkeeping below and try the next event.
        continue;
      }
      sign[idx_event] = f_t0 < 0 ? 1 : -1;
      found[idx_event] = true;
      found_any = true;
    }

    // If we found at least one event, then reset the solver state
    // back to the point of the event and apply all the events in
    // turn.
    if (found_any) {
      internals.last.interpolate(t1, y_next_.data());
      for (size_t idx_event = 0; idx_event < events.size(); ++idx_event) {
        if (found[idx_event]) {
          events[idx_event].action(t1, sign[idx_event], y_next_.data());
          internals.events.push_back({t1, idx_event, sign[idx_event]});
        }
      }
      internals.last.t1 = t1;
    }

    return {found_any, t1};
  }

  void apply_zero_every(real_type t, real_type* y,
                        const zero_every_type<real_type>& zero_every,
                        ode::internals<real_type>& internals) {
    // internals is not const because we write to .last.c4
    if (zero_every.empty() || internals.last.h == 0) {
      return;
    }
    for (const auto& el : zero_every) {
      const auto period = el.first;
      const auto t_last_step = internals.last.t0;
      const int n = std::floor(t / period);
      const int m = std::floor(t_last_step / period);
      if (n > m) {
        const auto t_reset = n * period;
        if (t_reset == t) {
          for (const auto j : el.second) {
            y[j] = 0;
          }
        } else {
          auto tmp = internals.last.c4.begin();
          const auto& index = el.second;
          internals.last.interpolate(t_reset, index, tmp);
          for (size_t i = 0; i < el.second.size(); ++i) {
            y[index[i]] -= tmp[i];
          }
        }
      }
    }
  }

  // action at the final step is slightly different:
  void apply_zero_every_final(real_type t, real_type* y,
                              const zero_every_type<real_type>& zero_every,
                              ode::internals<real_type>& internals) {
    if (zero_every.empty() || internals.last.h == 0) {
      return;
    }
    for (const auto& el : zero_every) {
      const auto period = el.first;
      if (internals.last.h > period) {
        const auto t_last_step = internals.last.t0;
        const int n = std::ceil(t / period);
        const int m = std::ceil(t_last_step / period);
        if (n > m) {
          const auto t_reset = (n - 1) * period;

          auto tmp = internals.last.c4.begin();
          const auto& index = el.second;
          internals.last.interpolate(t_reset, index, tmp);
          for (size_t i = 0; i < el.second.size(); ++i) {
            y[index[i]] -= tmp[i];
          }
        }
      }
    }
  }

private:
  size_t n_variables_;
  size_t n_special_;
  std::vector<real_type> y_next_;
  std::vector<real_type> y_stiff_;
  std::vector<real_type> k2_;
  std::vector<real_type> k3_;
  std::vector<real_type> k4_;
  std::vector<real_type> k5_;
  std::vector<real_type> k6_;
  real_type facc1_;
  real_type facc2_;

  // Internal private constants that define the integration tableau
  static constexpr real_type C2 = 0.2;
  static constexpr real_type C3 = 0.3;
  static constexpr real_type C4 = 0.8;
  static constexpr real_type C5 = 8.0 / 9.0;
  static constexpr real_type A21 = 0.2;
  static constexpr real_type A31 = 3.0 / 40.0;
  static constexpr real_type A32 = 9.0 / 40.0;
  static constexpr real_type A41 = 44.0 / 45.0;
  static constexpr real_type A42 = -56.0 / 15.0;
  static constexpr real_type A43 = 32.0 / 9.0;
  static constexpr real_type A51 = 19372.0 / 6561.0;
  static constexpr real_type A52 = -25360.0 / 2187.0;
  static constexpr real_type A53 = 64448.0 / 6561.0;
  static constexpr real_type A54 = -212.0 / 729.0;
  static constexpr real_type A61 = 9017.0 / 3168.0;
  static constexpr real_type A62 = -355.0 / 33.0;
  static constexpr real_type A63 = 46732.0 / 5247.0;
  static constexpr real_type A64 = 49.0 / 176.0;
  static constexpr real_type A65 = -5103.0 / 18656.0;
  static constexpr real_type A71 = 35.0 / 384.0;
  static constexpr real_type A73 = 500.0 / 1113.0;
  static constexpr real_type A74 = 125.0 / 192.0;
  static constexpr real_type A75 = -2187.0 / 6784.0;
  static constexpr real_type A76 = 11.0 / 84.0;
  static constexpr real_type E1 = 71.0 / 57600.0;
  static constexpr real_type E3 = -71.0 / 16695.0;
  static constexpr real_type E4 = 71.0 / 1920.0;
  static constexpr real_type E5 = -17253.0 / 339200.0;
  static constexpr real_type E6 = 22.0 / 525.0;
  static constexpr real_type E7 = -1.0 / 40.0;
  // ---- DENSE OUTPUT OF SHAMPINE (1986)
  static constexpr real_type D1 = -12715105075.0 / 11282082432.0;
  static constexpr real_type D3 = 87487479700.0 / 32700410799.0;
  static constexpr real_type D4 = -10690763975.0 / 1880347072.0;
  static constexpr real_type D5 = 701980252875.0 / 199316789632.0;
  static constexpr real_type D6 = -1453857185.0 / 822651844.0;
  static constexpr real_type D7 = 69997945.0 / 29380423.0;
};

}
}
