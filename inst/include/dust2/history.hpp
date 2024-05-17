#pragma once

#include <algorithm>

namespace dust {

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
      throw std::runtime_error("reorder extraction not implemented");
      /*
      // Default index:
      std::vector<size_t> index(n_particles_);
      for (size_t i = 0, k = 0; i < n_groups_; ++i) {
        for (size_t j = 0; j < n_particles_; ++j, ++k) {
          index_particle[k] = j;
        }
      }
      // the bits of this loop are independent, so just sort out the
      // bits to do a single particle at once or a time at once and
      // that will make it much easier to follow the core algorithm.
      for (size_t k = 0; k < position_; ++k) {
        const auto i = position_ - k - 1;
        const auto iter_order = order_.begin() + i * len_order_;
        const auto iter_state = state_.begin() + i * len_state_;
        const auto reorder = reorder[i];
        for (size_t ii = 0, kk = 0; ii < n_groups_; ++ii) {
          const auto offset_state = ii * n_state_ * n_particles_;
          const auto offset_order = ii * n_particles_;
          for (size_t jj = 0; jj < n_particles_; ++jj, ++kk) {
            std::copy_n(it_state + offset_state + index_particle[jj] * n_state_,
                        n_state_,
                        iter + offset_state + jj * n_state_);
            if (reorder) {
              index_particle[jj] =
                *(it_order + offset_order + index_particle[[jj]]);
            }
          }
        }
      }
      */
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
  std::vector<real_type> order_;
  std::vector<bool> reorder_;
};

}
