#pragma once

#include <deque>
#include <vector>
#include <dust2/tools.hpp>

namespace dust2 {
namespace ode {

// A single piece of history
template <typename real_type>
struct history_step {
  real_type t;
  real_type h;
  std::vector<real_type> c1;
  std::vector<real_type> c2;
  std::vector<real_type> c3;
  std::vector<real_type> c4;
  std::vector<real_type> c5;

  history_step(size_t n_variables) :
    c1(n_variables),
    c2(n_variables),
    c3(n_variables),
    c4(n_variables),
    c5(n_variables) {
  }

  history_step() {
  }

  history_step(real_type t, real_type h, std::vector<real_type> c1, std::vector<real_type> c2, std::vector<real_type> c3, std::vector<real_type> c4, std::vector<real_type> c5) :
    t(t),
    h(h),
    c1(c1),
    c2(c2),
    c3(c3),
    c4(c4),
    c5(c5) {
  }

  history_step subset(std::vector<size_t> index) const {
    return history_step(t,
                        h,
                        tools::subset(c1, index),
                        tools::subset(c2, index),
                        tools::subset(c3, index),
                        tools::subset(c4, index),
                        tools::subset(c5, index));
  }
};

template <typename real_type>
class history {
public:
  history(size_t n_state) :
    n_state_(n_state),
    n_state_total_(n_state_),
    initial_state_(n_state_) {
  }

  bool empty() const {
    return data_.empty();
  }

  size_t size() const {
    return data_.size();
  }

  size_t n_state() const {
    return n_state_;
  }

  // void set_index(const std::vector<size_t>& index) {
  //   if (tools::is_trivial_index(index, n_state_total_)) {
  //     index_.clear();
  //     n_state_ = n_state_total_;
  //   } else {
  //     index_ = index;
  //     n_state_ = index.size();
  //   }
  //   initial_state_.resize(n_state_);
  // }

  template <typename Iter>
  void reset(Iter state) {
    if (index_.empty()) {
      std::copy_n(state, n_state_, initial_state_.begin());
    } else {
      for (size_t i = 0; i < index_.size(); ++i) {
        initial_state_[i] = state[index_[i]];
      }
    }
    data_.clear();
  }

  void add(const history_step<real_type>& h) {
    if (index_.empty()) {
      data_.push_back(h);
    } else {
      data_.push_back(h.subset(index_));
    }
  }

  template <typename Iter>
  bool interpolate(real_type time, const std::vector<size_t>& index,
                   Iter result) {
    if (data_.empty() || time < data_[0].t) {
      for (size_t i = 0; i < index.size(); ++i) {
        *result = initial_state_[i];
        ++result;
      }
      return false;
    }

    const auto it = std::lower_bound(data_.begin(), data_.end(), time,
                                     [](auto value, auto& h) { value < h.t; });
    const auto u = (time - it->t) / it->h;
    const auto v = 1 - u;
    for (size_t i = 0; i < index.size(); ++i) {
      *result =
        it->c1[i] + u * (it->c2[i] + v * (it->c3[i] + u * (it->c4[i] + v * it->c5[i])));
      ++result;
    }
    return true;
  }

  auto& data() const {
    return data_;
  }

private:
  size_t n_state_;
  size_t n_state_total_;
  std::vector<size_t> index_;
  std::vector<real_type> initial_state_;
  std::deque<history_step<real_type>> data_;
};

}
}
