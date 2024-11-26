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

template <typename real_type>
struct event {
  using test_type = std::function<real_type(const real_type, const real_type*)>;
  using action_type = std::function<void(const real_type, const real_type, real_type*)>;
  std::vector<size_t> index;
  root_type root = root_type::both;
  test_type test;
  action_type action;

  event(const std::vector<size_t>& index, test_type test, action_type action, root_type root = root_type::both) :
    index(index), root(root), test(test), action(action) {
  }

  event(size_t index, action_type action, root_type root = root_type::both) :
    event({index}, [](real_type t, const real_type* y) { return y[0]; }, action, root) {
  }
};

template <typename real_type>
using events_type = std::vector<event<real_type>>;

template <typename real_type>
struct event_history_element {
  real_type time;
  size_t index;
  real_type sign;
};

template <typename real_type>
using event_history = std::vector<event_history_element<real_type>>;

}
}
