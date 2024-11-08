#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <dust2/tools.hpp>
#include <dust2/zero.hpp>
#include <dust2/continuous/control.hpp>
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

  std::vector<real_type> dydt;
  std::vector<real_type> step_times;
  real_type step_size;
  real_type error;
  size_t n_steps;
  size_t n_steps_accepted;
  size_t n_steps_rejected;

  internals(size_t n_variables) :
    last(n_variables),
    history_values(n_variables),
    dydt(n_variables) {
    reset();
  }

  void save_history() {
    history_values.add(last);
  }

  void reset() {
    step_size = 0;
    error = 0;
    n_steps = 0;
    n_steps_accepted = 0;
    n_steps_rejected = 0;
    step_times.clear();
    history_values.clear();
  }
};

template <typename real_type>
class solver {
public:
  solver(size_t n_variables, ode::control<real_type> control) :
    n_variables_(n_variables),
    control_(control),
    y_next_(n_variables_),
    y_stiff_(n_variables_),
    k2_(n_variables_),
    k3_(n_variables_),
    k4_(n_variables_),
    k5_(n_variables_),
    k6_(n_variables_),
    facc1_(1 / control_.factor_min),
    facc2_(1 / control_.factor_max) {
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
    const auto atol = control_.atol;
    const auto rtol = control_.rtol;
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
                 ode::internals<real_type>& internals, Rhs rhs) {
    auto success = false;
    auto reject = false;
    auto truncated = false;
    auto h = internals.step_size;

    while (!success) {
      if (internals.n_steps > control_.max_steps) {
        // throw a nicer error for all of these, with the current
        // time etc.
        throw std::runtime_error("too many steps");
      }
      if (h < control_.step_size_min) {
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
      const auto fac11 = std::pow(err, control_.constant);

      if (err <= 1) {
        success = true;
        accept(t, h, y, internals);
        internals.n_steps_accepted++;
        if (control_.debug_record_step_times) {
          internals.step_times.push_back(truncated ? t_end : t + h);
        }
        if (control_.save_history) {
          internals.save_history();
        }
        if (!truncated) {
          const auto fac_old =
            std::max(internals.error, static_cast<real_type>(1e-4));
          auto fac = fac11 / std::pow(fac_old, control_.beta);
          fac = clamp(fac / control_.factor_safe, facc2_, facc1_);
          const auto h_new = h / fac;
          const auto h_max = reject ? h : control_.step_size_max;
          internals.step_size = std::min(h_new, h_max);
          internals.error = err;
        }
      } else {
        reject = true;
        if (internals.n_steps_accepted >= 1) {
          internals.n_steps_rejected++;
        }
        h /= std::min(facc1_, fac11 / control_.factor_safe);
      }
    }

    return truncated ? t_end : t + h;
  }

  template <typename Rhs>
  void run(real_type t, real_type t_end, real_type* y,
           zero_every_type<real_type>& zero_every,
           ode::internals<real_type>& internals, Rhs rhs) {
    if (control_.critical_times.empty()) {
      while (t < t_end) {
        apply_zero_every(t, y, zero_every, internals);
        t = step(t, t_end, y, internals, rhs);
      }
    } else {
      // Slightly more complex loop which ensures we never integrate
      // over the times within our critical times.  The upper loop is
      // a special case of this but is kept simple.
      auto tc_end = control_.critical_times.end();
      auto tc = std::upper_bound(control_.critical_times.begin(), tc_end, t);
      auto t_end_i = (tc == tc_end || *tc >= t_end) ? t_end : *tc;
      while (t < t_end) {
        apply_zero_every(t, y, zero_every, internals);
        t = step(t, t_end_i, y, internals, rhs);
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
    internals.reset();
    if (control_.debug_record_step_times) {
      internals.step_times.push_back(t);
    }
    auto f0 = internals.dydt.data();
    auto f1 = k3_.data();
    auto y1 = y_next_.data();

    // Compute a first guess for explicit Euler as
    //   h = 0.01 * norm (y0) / norm (f0)
    // the increment for explicit euler is small compared to the solution
    rhs(t, y, f0);

    real_type norm_f = 0.0;
    real_type norm_y = 0.0;

    for (size_t i = 0; i < n_variables_; ++i) {
      const real_type sk = control_.atol + control_.rtol * std::abs(y[i]);
      norm_f += square(f0[i] / sk);
      norm_y += square(y[i] / sk);
    }
    // Magic numbers here, from Hairer
    real_type h = (norm_f <= 1e-10 || norm_y <= 1e-10) ?
      1e-6 : std::sqrt(norm_y / norm_f) * 0.01;
    h = std::min(h, control_.step_size_max);

    // Perform an explicit Euler step
    for (size_t i = 0; i < n_variables_; ++i) {
      y1[i] = y[i] + h * f0[i];
    }
    rhs(t + h, y1, f1);

    // Estimate the second derivative of the solution:
    real_type der2 = 0.0;
    for (size_t i = 0; i < n_variables_; ++i) {
      const real_type sk = control_.atol + control_.rtol * std::abs(y[i]);
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
    h = std::max(h, control_.step_size_min * 100);
    h = std::min(h, control_.step_size_max);
    internals.step_size = h;
  }

private:
  void accept(real_type t, real_type h, real_type* y, ode::internals<real_type>& internals) {
    // We might want to only do this bit if we'll actually use the
    // history, but it's pretty cheap really.
    internals.last.t = t;
    internals.last.h = h;
    for (size_t i = 0; i < n_variables_; ++i) {
      const auto ydiff = y_next_[i] - y[i];
      const auto bspl = h * internals.dydt[i] - ydiff;
      internals.last.c1[i] = y[i];
      internals.last.c2[i] = ydiff;
      internals.last.c3[i] = bspl;
      internals.last.c4[i] = -h * k2_[i] + ydiff - bspl;
    }

    std::copy_n(k2_.begin(), n_variables_, internals.dydt.begin());
    std::copy_n(y_next_.begin(), n_variables_, y);
  }

  real_type interpolate(size_t idx,
                        const real_type theta,
                        const real_type theta1,
                        const ode::internals<real_type>& internals) {
    return internals.last.c1[idx] + theta *
      (internals.last.c2[idx] + theta1 *
       (internals.last.c3[idx] + theta *
        (internals.last.c4[idx] + theta1 *
         internals.last.c5[idx])));
  }

  void apply_zero_every(real_type t, real_type* y,
                        const zero_every_type<real_type>& zero_every,
                        const ode::internals<real_type>& internals) {
    if (zero_every.empty() || internals.last.h == 0) {
      return;
    }
    for (const auto& el : zero_every) {
      const auto period = el.first;
      const auto t_last_step = t - internals.last.h;
      const int n = std::floor(t / period);
      const int m = std::floor(t_last_step / period);
      if (n > m) {
        const auto t_reset = n * period;
        if (t_reset == t) {
          for (const auto j : el.second) {
            y[j] = 0;
          }
        } else {
          const auto theta = (t_reset - t_last_step) / internals.last.h;
          const auto theta1 = 1 - theta;
          for (const auto j : el.second) {
            y[j] -= interpolate(j, theta, theta1, internals);
          }
        }
      }
    }
  }

  // action at the final step is slightly different:
  void apply_zero_every_final(real_type t, real_type* y,
                              const zero_every_type<real_type>& zero_every,
                              const ode::internals<real_type>& internals) {
    if (zero_every.empty() || internals.last.h == 0) {
      return;
    }
    for (const auto& el : zero_every) {
      const auto period = el.first;
      if (internals.last.h > period) {
        const auto t_last_step = t - internals.last.h;
        const int n = std::ceil(t / period);
        const int m = std::ceil(t_last_step / period);
        if (n > m) {
          const auto t_reset = (n - 1) * period;
          const auto theta = (t_reset - t_last_step) / internals.last.h;
          const auto theta1 = 1 - theta;
          for (const auto j : el.second) {
            y[j] -= interpolate(j, theta, theta1, internals);
          }
        }
      }
    }
  }

  size_t n_variables_;
  ode::control<real_type> control_;
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
