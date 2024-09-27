#pragma once

#include <dust2/r/helpers.hpp>
#include <dust2/continuous/control.hpp>

namespace dust2 {
namespace r {

template <typename real_type>
dust2::ode::control<real_type> validate_ode_control(cpp11::list r_time_control) {
  cpp11::list ode_control = r_time_control["ode_control"];
  const auto max_steps = dust2::r::read_int(ode_control, "max_steps");
  const auto atol = dust2::r::read_real(ode_control, "atol");
  const auto rtol = dust2::r::read_real(ode_control, "rtol");
  const auto step_size_min =
    std::max(static_cast<real_type>(dust2::r::read_real(ode_control, "step_size_min")),
             std::numeric_limits<real_type>::epsilon());
  const auto step_size_max = dust2::r::read_real(ode_control, "step_size_max");
  const auto debug_record_step_times =
    dust2::r::read_bool(ode_control, "debug_record_step_times");
  return dust2::ode::control<real_type>(max_steps, atol, rtol, step_size_min,
                                        step_size_max, debug_record_step_times);
}

}
}
