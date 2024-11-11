#pragma once

#include <algorithm>
#include <set>
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
    // std::set<size_t> tmp;
    // for (auto& el : delays_) {
    //   tmp.insert(el.index.begin(), el.index.end());
    // }
    // index_state_.insert(tmp.begin(), tmp.end());
    // std::sort(index_state_.begin(), index_state_.end());
    // for (size_t i = 0; i < delays_.size(); ++i) {
    //   for (auto j : delays_[i].index) {
    //     const auto k = *std::find(index_state_.begin(), index_state_.end(), j);
    //     index_out_[i].push_back(k);
    //   }
    // }
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
  std::vector<size_t> index_state_;
  std::vector<std::vector<size_t>> index_out_;
};


}
}
