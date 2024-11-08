#pragma once

namespace dust2 {
namespace ode{

// A single piece of history
template <typename real_type>
struct history {
  real_type t;
  real_type h;
  std::vector<real_type> c1;
  std::vector<real_type> c2;
  std::vector<real_type> c3;
  std::vector<real_type> c4;
  std::vector<real_type> c5;

  history(size_t n_variables) :
    c1(n_variables),
    c2(n_variables),
    c3(n_variables),
    c4(n_variables),
    c5(n_variables) {
  }

  history() {
  }

  history(real_type t, real_type h, std::vector<real_type> c1, std::vector<real_type> c2, std::vector<real_type> c3, std::vector<real_type> c4, std::vector<real_type> c5) :
    t(t),
    h(h),
    c1(c1),
    c2(c2),
    c3(c3),
    c4(c4),
    c5(c5) {
  }

  history subset(std::vector<size_t> index) {
    return history(t,
                   h,
                   tools::subset(c1, index),
                   tools::subset(c2, index),
                   tools::subset(c3, index),
                   tools::subset(c4, index),
                   tools::subset(c5, index));
  }
};

}
}
