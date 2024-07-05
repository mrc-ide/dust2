#pragma once

#include <map>
#include <vector>

namespace dust2 {

template <typename real_type>
using zero_every_type = std::map<real_type, std::vector<size_t>>;

template <typename T>
auto zero_every(const typename T::shared_state& shared) {
  using real_type = typename T::real_type;
  return zero_every_type<real_type>();
}

template <typename T>
auto zero_every_vec(const std::vector<typename T::shared_state>& shared) {
  using real_type = typename T::real_type;
  std::vector<zero_every_type<real_type>> ret;
  ret.reserve(shared.size());
  for (const auto& el : shared) {
    ret.push_back(zero_every<T>(el));
  }
  return ret;
}

}
