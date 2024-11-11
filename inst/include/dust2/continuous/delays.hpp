#pragma once

#include <vector>

#include <dust2/continuous/history.hpp>

namespace dust2 {
namespace ode {

template <typename real_type>
struct delay {
  real_type tau;
  std::vector<size_t> index;
};

template <typename real_type>
using delay_result_type = std::vector<std::vector<real_type>>;

template <typename real_type>
class delays {
public:
  delays(std::initializer_list<delay<real_type>> data) : delays_(data) {
  }

  delay_result_type<real_type> result() {
    delay_result_type<real_type> ret;
    ret.reserve(size());
    for (const auto& el : delays_) {
      ret.push_back(std::vector<real_type>(el.index.size()));
    }
    return ret;
  }

  void eval(real_type time, const history<real_type>& h,
            delay_result_type<real_type>& result) const {
    for (size_t i = 0; i < delays_.size(); ++i) {
      h.interpolate(time - delays_[i].tau,
                    delays_[i].index,
                    result[i].begin());
    }
  }

  size_t size() const {
    return delays_.size();
  }

private:
  std::vector<delay<real_type>> delays_;
};


}
}
