#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <numeric>

namespace dust2 {
namespace tools {

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

template <typename T, typename U>
T accumulate_periodic(T time, T period, U previous, U value) {
  return is_evenly_divisible_by(time, period) ? value : previous + value;
}

inline size_t thread_index() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

template <typename T>
bool all(const T& x) {
  return std::all_of(x.begin(), x.end(), [](auto v) { return v; });
}

template <typename T>
bool any(const T& x) {
  return std::any_of(x.begin(), x.end(), [](auto v) { return v; });
}

template <typename T>
T prod(const std::vector<T>& x) {
  return std::accumulate(x.begin(), x.end(), 1, std::multiplies<>{});
}

}
}
