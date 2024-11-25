#pragma once

#include <functional>
#include <vector>

namespace dust2 {
namespace ode {

// Do we detect an event when passing through the root as we increase
// (negative to positive), decrease (positive to negative) or both:
enum class root_type {
  both,
  increase,
  decrease
};

// The actual logic for the above
template <typename real_type>
bool is_root(const real_type a, const real_type b, const root_type& root) {
  switch(root) {
  case root_type::both:
    return a * b < 0;
  case root_type::increase:
    return a < 0 && b > 0;
  case root_type::decrease:
    return a > 0 && b < 0;
  }
  return false;
}

template<typename real_type>
struct event {
  size_t index;
  real_type value;
  std::function<void(real_type, real_type*)> action; // time, y
  root_type root = root_type::both;
};

template<typename real_type>
using events_type = std::vector<event<real_type>>;

}
}
