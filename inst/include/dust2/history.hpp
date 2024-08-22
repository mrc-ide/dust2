#pragma once

#include <algorithm>
#include <stdexcept>
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
  void add(real_type time, IterReal iter_state,
           const std::vector<size_t>& index_group) {
    copy_state_(iter_state, index_group);
    update_position(time, false);
  }

  template <typename IterReal, typename IterSize>
  void add(real_type time, IterReal iter_state, IterSize iter_order,
           const std::vector<size_t>& index_group) {
    copy_state_(iter_state, index_group);
    copy_order_(iter_order, index_group);
    update_position(time, true);
  }

  template <typename IterReal, typename IterSize>
  void add_with_index(real_type time, IterReal iter_state, IterSize iter_index,
                      size_t n_state_total,
                      const std::vector<size_t>& index_group) {
    copy_state_with_index_(iter_state, iter_index, n_state_total, index_group);
    update_position(time, false);
  }

  template <typename IterReal, typename IterSize>
  void add_with_index(real_type time, IterReal iter_state, IterSize iter_order,
                      IterSize iter_index, size_t n_state_total,
                      const std::vector<size_t>& index_group) {
    copy_state_with_index_(iter_state, iter_index, n_state_total, index_group);
    copy_order_(iter_order, index_group);
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

  auto size_order() const {
    return position_ * len_order_;
  }

  auto n_state() const {
    return n_state_;
  }

  auto n_times() const {
    return position_;
  }

  template <typename Iter>
  void export_time(Iter iter) const {
    std::copy_n(times_.begin(), position_, iter);
  }

  template <typename Iter>
  void export_state(Iter iter, bool reorder,
                    const std::vector<size_t>& index_group,
		    const std::vector<size_t>& select_particle) const {
    reorder = reorder && n_particles_ > 1 && position_ > 0 &&
      // tools::any(reorder_);
      std::any_of(reorder_.begin(), reorder_.end(), [](auto v) { return v; });

    auto use_select_particle = !select_particle.empty();
    auto use_index_group = index_group.size() < n_groups_;
    if (!use_index_group) {
      for (size_t i = 0; i < n_groups_; ++i) {
        if (index_group[i] != i) {
          use_index_group = true;
          break;
        }
      }
    }

    if (!reorder && !use_index_group && !use_select_particle) {
      std::copy_n(state_.begin(), position_ * len_state_, iter);
    } else if (reorder) {
      const size_t n_particles_out = use_select_particle ? 1 : n_particles_;
      // Default index - we can improve this by fixing k to be an
      // offset, and also by saving this within this history object.
      std::vector<size_t> index_particle(n_particles_out * n_groups_);
      for (size_t i = 0, k = 0; i < n_groups_; ++i) {
	if (use_select_particle) {
	  index_particle[i] = select_particle[i];
	} else {
	  for (size_t j = 0; j < n_particles_; ++j, ++k) {
	    index_particle[k] = j;
	  }
	}
      }

      const auto len_state_output =
        n_state_ * n_particles_out * index_group.size();
      for (size_t irev = 0; irev < position_; ++irev) {
        const auto i = position_ - irev - 1;
        const auto iter_order = order_.begin() + i * len_order_;
        const auto iter_state = state_.begin() + i * len_state_;
        for (size_t j = 0; j < index_group.size(); ++j) {
          const auto k = index_group[j];
          const auto offset_state = k * n_state_ * n_particles_;
          const auto offset_index = k * n_particles_out;
          const auto offset_output =
            i * len_state_output + j * n_state_ * n_particles_out;
	  if (use_select_particle) {
	    reorder_single_(iter_state + offset_state,
			    iter_order + offset_index,
			    reorder_[i],
			    iter + offset_output,
			    index_particle[j]);
	  } else {
	    reorder_group_(iter_state + offset_state,
			   iter_order + offset_index,
			   reorder_[i],
			   iter + offset_output,
			   index_particle.begin() + offset_index);
	  }
        }
      }
    } else {
      // use_index_group || use_select_particle (or both!)
      for (size_t i = 0; i < position_; ++i) {
	const auto iter_state = state_.begin() + i * len_state_;
	for (auto j : index_group) {
	  const auto offset_state = j * n_state_ * n_particles_;
	  if (use_select_particle) {
	    const auto offset_j = select_particle[j] * n_state_;
	    iter = std::copy_n(iter_state + offset_state + offset_j,
			       n_state_,
			       iter);
	  } else {
	    iter = std::copy_n(iter_state + offset_state,
			       n_state_ * n_particles_,
			       iter);
	  }
	}
      }
    }
  }

  template <typename Iter>
  void export_order(Iter iter) const {
    std::copy_n(order_.begin(), position_ * len_order_, iter);
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

  // Reference implementation for this is mcstate::history_single and
  // mcstate::history_multiple
  template <typename Iter>
  void reorder_group_(typename std::vector<real_type>::const_iterator iter_state,
                      typename std::vector<size_t>::const_iterator iter_order,
                      bool reorder,
                      Iter iter_dest,
                      typename std::vector<size_t>::iterator index_particle) const {
    for (size_t i = 0; i < n_particles_; ++i) {
      std::copy_n(iter_state + *(index_particle + i) * n_state_,
                  n_state_,
                  iter_dest + i * n_state_);
      if (reorder) {
        const auto index_particle_i = index_particle + i;
        *index_particle_i = *(iter_order + *index_particle_i);
      }
    }
  }

  template <typename Iter>
  void reorder_single_(typename std::vector<real_type>::const_iterator iter_state,
		       typename std::vector<size_t>::const_iterator iter_order,
		       bool reorder,
		       Iter iter_dest,
		       size_t& index_particle) const {
    std::copy_n(iter_state + index_particle * n_state_,
		n_state_,
		iter_dest);
    if (reorder) {
      index_particle = *(iter_order + index_particle);
    }
  }

  void update_position(real_type time, bool reorder) {
    times_[position_] = time;
    reorder_[position_] = reorder;
    position_++;
  }

  template <typename IterReal>
  void copy_state_(IterReal iter, const std::vector<size_t>& index_state) {
    std::copy_n(iter, len_state_, state_.begin() + len_state_ * position_);
  }

  template <typename IterSize>
  void copy_order_(IterSize iter_order,
                   const std::vector<size_t>& index_state) {
    std::copy_n(iter_order, len_order_,
                order_.begin() + position_ * len_order_);
  }

  template <typename IterReal, typename IterSize>
  void copy_state_with_index_(IterReal iter_state,
                              IterSize iter_index_state,
                              size_t n_state_total,
                              const std::vector<size_t>& index_group) {
    auto iter_dest = state_.begin() + position_ * len_state_;
    for (size_t i = 0; i < n_groups_ * n_particles_; ++i) {
      const auto iter_state_i = iter_state + i * n_state_total;
      auto iter_index_state_i = iter_index_state;
      for (size_t k = 0; k < n_state_; ++k, ++iter_dest, ++iter_index_state_i) {
        *iter_dest = *(iter_state_i + *iter_index_state_i);
      }
    }
  }
};

}
