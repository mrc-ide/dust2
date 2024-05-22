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
    reorder_(n_times) {
  }

  void resize(size_t n_times) {
    n_times_ = n_times;
    state_.resize(len_state_ * n_times_);
    order_.resize(len_order_ * n_times_);
    reset();
  }

  void reset() {
    position_ = 0;
  }

  template <typename Iter>
  void add(real_type time, Iter iter) {
    // TODO: bounds check here (and in the method below)?
    std::copy_n(iter, len_state_, state_.begin() + len_state_ * position_);
    times_[position_] = time;
    reorder_[position_] = false;
    position_++;
  }

  template <typename IterReal, typename IterSize>
  void add(real_type time, IterReal iter_state, IterSize iter_order) {
    // This can't easily call add(Iter, real_type) because we need
    // read position_ and write reorder_; the duplication is minimal
    // though.
    std::copy_n(iter_state, len_state_,
                state_.begin() + position_ * len_state_);
    std::copy_n(iter_order, len_order_,
                order_.begin() + position_ * len_order_);
    times_[position_] = time;
    reorder_[position_] = true;
    position_++;
  }

  // These allow a consumer to allocate the right size structures for
  // time and state for the total that we've actually used.
  auto size_time() const {
    return position_;
  }

  auto size_state() const {
    return position_ * len_state_;
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
};

}
