#pragma once

#include <algorithm>
#include <set>
#include <vector>

#include <dust2/continuous/history.hpp>

namespace dust2 {
namespace ode {

struct delay_variable {
  size_t offset;
  size_t length;
};

template <typename real_type>
struct delay {
  real_type tau;
  std::vector<delay_variable> variables;
  // This is the index into the main state
  std::vector<size_t> index;
  // This is the offset to find the variable we are interested in,
  // within the data returned from the delays.
  std::vector<size_t> offset;
  delay(real_type tau, std::initializer_list<delay_variable> variables) :
    tau(tau), variables(variables) {
    size_t len = 0;
    for (auto& el : variables) {
      offset.push_back(len);
      for (size_t i = 0; i < el.length; ++i) {
        index.push_back(el.offset + i);
      }
      len += el.length;
    }
  }
};

template <typename real_type>
struct delay_result_element {
  std::vector<real_type> data;
  std::vector<size_t> offset;

  delay_result_element(size_t size, std::vector<size_t> offset) :
    data(size), offset(offset) {
  }
};

template <typename real_type>
using delay_result_type = std::vector<delay_result_element<real_type>>;

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

  auto result() const {
    delay_result_type<real_type> ret;
    ret.reserve(size());
    for (const auto& el : delays_) {
      ret.push_back({el.index.size(), el.offset});
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
                    index_out_[i],
                    result[i].data.begin());
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
