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
  delays(std::initializer_list<delay<real_type>> data) :
    delays_(data),
    index_out_(delays_.size()) {
    std::set<size_t> tmp;
    for (auto& el : delays_) {
      tmp.insert(el.index.begin(), el.index.end());
    }
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(index_save_));
    std::sort(index_save_.begin(), index_save_.end());
    for (size_t i = 0; i < delays_.size(); ++i) {
      for (auto j : delays_[i].index) {
        const auto it = std::find(index_save_.begin(), index_save_.end(), j);
        index_out_[i].push_back(std::distance(index_save_.begin(), it));
      }
    }
  }

  auto& index() const {
    return index_save_;
  }

  delay_result_type<real_type> result() {
    delay_result_type<real_type> ret;
    ret.reserve(size());
    for (const auto& el : delays_) {
      ret.push_back(std::vector<real_type>(el.index.size()));
    }
    return ret;
  }

  real_type step_size_max(real_type h, bool used_in_rhs) const {
    if (used_in_rhs) {
      h = std::accumulate(delays_.begin(), delays_.end(), h,
                          [](auto a, auto b) { return a < b.tau ? a : b.tau; });
    }
    return h;
  }

  void eval(real_type time, const history<real_type>& h,
            delay_result_type<real_type>& result) const {
    for (size_t i = 0; i < delays_.size(); ++i) {
      h.interpolate(time - delays_[i].tau,
                    index_out_[i], // or save back into delays_[i].index above?
                    result[i].begin());
    }
  }

  size_t size() const {
    return delays_.size();
  }

private:
  std::vector<delay<real_type>> delays_;
  std::vector<size_t> index_save_;
  std::vector<std::vector<size_t>> index_out_;
};


}
}
