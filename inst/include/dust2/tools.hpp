#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

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

inline std::vector<size_t> integer_sequence(size_t n) {
  std::vector<size_t> ret;
  ret.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    ret.push_back(i);
  }
  return ret;
}

template <typename T>
std::vector<T> subset(const std::vector<T>& x, const std::vector<size_t> index) {
  std::vector<T> ret;
  ret.reserve(index.size());
  for (auto i : index) {
    ret.push_back(x[i]);
  }
  return ret;
}

inline bool is_trivial_index(const std::vector<size_t>& index, size_t n) {
  if (index.empty()) {
    return true;
  }
  if (index.size() != n) {
    return false;
  }
  for (size_t i = 0; i < n; ++i) {
    if (index[i] != i) {
      return false;
    }
  }
  return true;
}

}
}
