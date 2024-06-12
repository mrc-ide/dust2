#pragma once

#include <cmath>

namespace dust2 {
namespace tools {

namespace {
template <typename T>
bool is_evenly_divisible_by(T num, T by);

template <>
inline bool is_evenly_divisible_by(double num, double by) {
  // This eps is chosen to be less than machine precision, though more
  // than sqrt(precision) and is consistent with the maximum expected
  // accumulation of rounding error for pathological choices of dt.
  constexpr double eps = 1e-13;
  return std::abs(std::fmod(num, by)) < eps;
}

}

template <typename T, typename U>
T accumulate_periodic(T time, T period, U previous, U value) {
  return is_evenly_divisible_by(time, period) ? value : previous + value;
}

}
}
