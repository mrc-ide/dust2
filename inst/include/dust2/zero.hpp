#pragma once

#include <map>
#include <vector>

namespace dust2 {

template <typename T>
using zero_every_type = std::map<typename T::real_type, std::vector<size_t>>;

template <typename T>
auto zero_every(const typename T::shared_state& shared) {
  return zero_every_type<T>();
}

template <typename T>
auto zero_every_vec(const std::vector<typename T::shared_state>& shared) {
  std::vector<zero_every_type<T>> ret;
  ret.reserve(shared.size());
  for (const auto& el : shared) {
    ret.push_back(zero_every<T>(el));
  }
  return ret;
}

}
