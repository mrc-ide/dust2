#pragma once

#include <algorithm>
#include <vector>

namespace dust2 {

// We might want a version of this that saves a subset of state too,
// we can think about that later though.  We also need a version that
// can allow working back through the graph of history.
template <typename real_type>
class history {
public:
  history(size_t n_state, size_t n_particles, size_t n_groups, size_t n_times) :
    n_state_(n_state),
    n_particles_(n_particles),
    n_groups_(n_groups),
    n_times_(n_times),
    len_state_(n_state_ * n_particles_ * n_groups_),
    len_order_(n_particles * n_groups_),
    position_(0),
    times_(n_times_),
    state_(len_state_ * n_times_),
    order_(len_order_ * n_times_),
    reorder_(n_times),
    dims_({n_state_, n_particles_, n_groups_, position_}) {
  }

  void resize_state(size_t n_state) {
    if (n_state_ != n_state) {
      n_state_ = n_state;
      len_state_ = (n_state_ * n_particles_ * n_groups_);
      state_.resize(len_state_ * n_times_);
    }
    reset();
  }

  void resize_time(size_t n_times) {
    if (n_times_ != n_times) {
      n_times_ = n_times;
      state_.resize(len_state_ * n_times_);
      order_.resize(len_order_ * n_times_);
    }
    reset();
  }

  void reset() {
    position_ = 0;
  }

  template <typename IterReal>
  void add(real_type time, IterReal iter_state) {
    copy_state_(iter_state);
    update_position(time, false);
  }

  template <typename IterReal, typename IterSize>
  void add(real_type time, IterReal iter_state, IterSize iter_order) {
    copy_state_(iter_state);
    copy_order_(iter_order);
    update_position(time, true);
  }

  template <typename IterReal, typename IterSize>
  void add_with_index(real_type time, IterReal iter_state, IterSize iter_index,
                      size_t n_state_total) {
    copy_state_with_index_(iter_state, iter_index, n_state_total);
    update_position(time, false);
  }

  template <typename IterReal, typename IterSize>
  void add_with_index(real_type time, IterReal iter_state, IterSize iter_order,
                      IterSize iter_index, size_t n_state_total) {
    copy_state_with_index(iter_state, iter_index, n_state_total);
    copy_order_(iter_order);
    update_position(time, true);
  }

  // These allow a consumer to allocate the right size structures for
  // time and state for the total that we've actually used.
  auto size_time() const {
    return position_;
  }

  auto size_state() const {
    return position_ * len_state_;
  }

  auto& dims() const {
    return dims_;
  }

  template <typename Iter>
  void export_time(Iter iter) {
    std::copy_n(times_.begin(), position_, iter);
  }

  template <typename Iter>
  void export_state(Iter iter, bool reorder) {
    reorder = reorder && n_particles_ > 1 && position_ > 0 &&
      std::any_of(reorder_.begin(), reorder_.end(), [](auto v) { return v; });
    if (reorder) {
      // Default index:
      std::vector<size_t> index_particle(n_particles_ * n_groups_);
      for (size_t i = 0, k = 0; i < n_groups_; ++i) {
        for (size_t j = 0; j < n_particles_; ++j, ++k) {
          index_particle[k] = j;
        }
      }

      for (size_t irev = 0; irev < position_; ++irev) {
        const auto i = position_ - irev - 1; // can move this to be the loop
        const auto iter_order = order_.begin() + i * len_order_;
        const auto iter_state = state_.begin() + i * len_state_;
        // This bit here is independent among groups
        for (size_t j = 0; j < n_groups_; j++) {
          const auto offset_state = j * n_state_ * n_particles_;
          const auto offset_index = j * n_particles_;
          reorder_group_(iter_state + offset_state,
                         iter_order + offset_index,
                         reorder_[i],
                         iter + i * len_state_ + offset_state,
                         index_particle.begin() + offset_index);
        }
      }
    } else {
      // No reordering is requested or possible so just dump out directly:
      std::copy_n(state_.begin(), position_ * len_state_, iter);
    }
  }

private:
  size_t n_state_;
  size_t n_particles_;
  size_t n_groups_;
  size_t n_times_;
  size_t len_state_; // length of an update to state
  size_t len_order_; // length of an update to order
  size_t position_;
  std::vector<real_type> times_;
  std::vector<real_type> state_;
  std::vector<size_t> order_;
  std::vector<bool> reorder_;
  std::array<size_t, 4> dims_;

  // Reference implementation for this is mcstate:::history_single and
  // mcstate::history_multiple
  template <typename Iter>
  void reorder_group_(typename std::vector<real_type>::const_iterator iter_state,
                      typename std::vector<size_t>::const_iterator iter_order,
                      bool reorder,
                      Iter iter_dest,
                      typename std::vector<size_t>::iterator index) {
    for (size_t i = 0; i < n_particles_; ++i) {
      std::copy_n(iter_state + *(index + i) * n_state_,
                  n_state_,
                  iter_dest + i * n_state_);
      if (reorder) {
        const auto index_i = index + i;
        *index_i = *(iter_order + *index_i);
      }
    }
  }

  void update_position(real_type time, bool reorder) {
    times_[position_] = time;
    reorder_[position_] = reorder;
    position_++;
    dims_[3] = position_;
  }

  template <typename IterReal>
  void copy_state_(IterReal iter) {
    std::copy_n(iter, len_state_, state_.begin() + len_state_ * position_);
  }

  template <typename IterSize>
  void copy_order_(IterSize iter_order) {
    std::copy_n(iter_order, len_order_,
                order_.begin() + position_ * len_order_);
  }

  template <typename IterReal, typename IterSize>
  void copy_state_with_index_(IterReal iter_state, IterSize iter_index,
                              size_t n_state_total) {
    auto iter_dest = state_.begin() + position_ * len_state_;
    for (size_t i = 0; i < n_groups_ * n_particles_; ++i) {
      const auto iter_state_i = iter_state + i * n_state_total;
      auto iter_index_i = iter_index;
      for (size_t k = 0; k < n_state_; ++k, ++iter_dest, ++iter_index_i) {
        *iter_dest = *(iter_state_i + *iter_index_i);
      }
    }
  }
};

}
