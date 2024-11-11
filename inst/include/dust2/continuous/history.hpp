#pragma once

#include <deque>
#include <stdexcept>
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
    t(0),
    h(0),
    c1(n_variables),
    c2(n_variables),
    c3(n_variables),
    c4(n_variables),
    c5(n_variables) {
  }

  history_step() : t(0), h(0) {
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

  template <typename Iter>
  void interpolate(real_type time, const std::vector<size_t>& index,
                   Iter result) const {
    const auto u = (time - t) / h;
    const auto v = 1 - u;
    for (auto i : index) {
      *result = c1[i] + u * (c2[i] + v * (c3[i] + u * (c4[i] + v * c5[i])));
      ++result;
    }
  }

  template <typename Iter>
  void interpolate(real_type time, Iter result) const {
    const auto u = (time - t) / h;
    const auto v = 1 - u;
    for (size_t i = 0; i < c1.size(); ++i) {
      *result = c1[i] + u * (c2[i] + v * (c3[i] + u * (c4[i] + v * c5[i])));
      ++result;
    }
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

  // We may enable this once everything is working, but most likely
  // it's just coming in via the constructor.  Here so I remember it.
  //
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
                   Iter result) const {
    if (data_.empty() || time < data_[0].t) {
      for (size_t i = 0; i < index.size(); ++i) {
        *result = initial_state_[i];
        ++result;
      }
      return false;
    } else {
      const auto it = std::lower_bound(data_.begin(), data_.end(), time,
                                       [](const auto& h, const auto& value) {
                                         return value > (h.t + h.h);
                                       });
      if (it == data_.end()) {
        throw std::runtime_error("Failed to locate history element");
      }
      it->interpolate(time, index, result);
      return true;
    }
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
