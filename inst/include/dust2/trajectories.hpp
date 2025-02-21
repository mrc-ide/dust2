#pragma once

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <dust2/tools.hpp>

namespace dust2 {

template <typename real_type>
class trajectories {
public:
  trajectories(size_t n_state, size_t n_particles, size_t n_groups, size_t n_times) :
    n_state_(n_state),
    n_state_total_(n_state),
    n_particles_(n_particles),
    n_groups_(n_groups),
    n_groups_total_(n_groups),
    n_times_(n_times),
    len_state_(n_state_total_ * n_particles_ * n_groups_total_),
    len_order_(n_particles_ * n_groups_total_),
    position_(0) {
  }

  void set_index_and_reset(const std::vector<size_t>& index_state,
                           const std::vector<size_t>& index_group) {
    index_state_ = index_state;
    index_group_ = index_group;
    use_index_state_ = !tools::is_trivial_index(index_state_, n_state_total_);
    use_index_group_ = !tools::is_trivial_index(index_group_, n_groups_total_);
    n_state_ = use_index_state_ ? index_state.size() : n_state_total_;
    n_groups_ = use_index_group_ ? index_group.size() : n_groups_total_;
    len_state_ = n_state_ * n_particles_ * n_groups_;

    times_.resize(n_times_);
    state_.resize(len_state_ * n_times_);
    order_.resize(len_order_ * n_times_);
    reorder_.resize(n_times_);

    position_ = 0;
  }

  void set_n_times_and_reset(size_t n_times,
                             const std::vector<size_t>& index_group) {
    n_times_ = n_times;
    set_index_and_reset(std::vector<size_t>(), index_group);
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

  template <typename IterSize>
  void update_order(IterSize iter_order) {
    // TODO - implementation here.
  }

  // These allow a consumer to allocate the right size structures for
  // time and state for the total that we've actually used.
  auto size_time() const {
    return position_;
  }

  auto n_state() const {
    return n_state_;
  }

  auto n_particles() const {
    return n_particles_;
  }

  auto n_groups() const {
    return n_groups_;
  }

  auto n_times() const {
    return position_;
  }

  auto& index_group() const {
    return index_group_;
  }

  template <typename Iter>
  void export_time(Iter iter) const {
    std::copy_n(times_.begin(), position_, iter);
  }

  // Note that this uses 'select_particle' to avoid conflict with our
  // internal 'index_particle'. Also we are not extracting the
  // rectangle defined by the combination of the group/particle values
  // but pairs of (group,particle)'s.
  template <typename Iter>
  void export_state(Iter iter, bool reorder,
                    const std::vector<size_t>& select_particle) const {
    reorder = reorder && n_particles_ > 1 && position_ > 0 &&
      tools::any(reorder_);

    auto use_select_particle = !select_particle.empty();

    if (!reorder && !use_select_particle) {
      // Optimised for simplest case of just dumping out everything
      std::copy_n(state_.begin(), position_ * len_state_, iter);
    } else if (reorder) {
      const size_t n_particles_out = use_select_particle ? 1 : n_particles_;
      std::vector<size_t> index_particle(n_particles_out * n_groups_);
      for (size_t i = 0; i < n_groups_; ++i) {
        if (use_select_particle) {
          index_particle[i] = select_particle[i];
        } else {
          for (size_t j = 0, k = i * n_particles_; j < n_particles_; ++j, ++k) {
            index_particle[k] = j;
          }
        }
      }

      const auto len_state_out = n_state_ * n_particles_out * n_groups_;
      for (size_t irev = 0; irev < position_; ++irev) {
        const auto i = position_ - irev - 1;
        const auto iter_order = order_.begin() + i * len_order_;
        const auto iter_state = state_.begin() + i * len_state_;
        for (size_t j = 0; j < n_groups_; ++j) {
          const auto offset_state = j * n_state_ * n_particles_;
          const auto offset_index = j * n_particles_;
          const auto offset_output =
            i * len_state_out + j * n_state_ * n_particles_out;
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
      for (size_t i = 0; i < position_; ++i) {
        const auto iter_state = state_.begin() + i * len_state_;
        for (size_t j = 0; j < n_groups_; ++j) {
          const auto offset_state = j * n_state_ * n_particles_;
          const auto offset_j = select_particle[j] * n_state_;
          iter = std::copy_n(iter_state + offset_state + offset_j,
                             n_state_,
                             iter);
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
  size_t n_state_total_;
  size_t n_particles_;
  size_t n_groups_;
  size_t n_groups_total_;
  size_t n_times_;
  size_t len_state_; // length of an update to state
  size_t len_order_; // length of an update to order
  size_t position_;
  std::vector<real_type> times_;
  std::vector<real_type> state_;
  std::vector<size_t> order_;
  std::vector<bool> reorder_;
  std::vector<size_t> index_state_;
  std::vector<size_t> index_group_;
  bool use_index_state_;
  bool use_index_group_;

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
  void copy_state_(IterReal iter_src) {
    auto iter_dst = state_.begin() + position_ * len_state_;
    if (!use_index_state_ && !use_index_group_) {
      std::copy_n(iter_src, len_state_, iter_dst);
    } else if (!use_index_state_) {
      const auto len = n_state_ * n_particles_;
      for (auto i : index_group_) {
        const auto iter_src_i = iter_src + i * n_particles_ * n_state_total_;
        iter_dst = std::copy_n(iter_src_i, len, iter_dst);
      }
    } else { // we always end up here if we have a state index
      for (size_t i = 0; i < n_groups_; ++i) {
        for (size_t j = 0; j < n_particles_; ++j) {
          const auto iter_src_i = iter_src +
            (index_group_[i] * n_particles_ + j) * n_state_total_;
          for (size_t k = 0; k < n_state_; ++k, ++iter_dst) {
            *iter_dst = *(iter_src_i + index_state_[k]);
          }
        }
      }
    }
  }

  template <typename IterSize>
  void copy_order_(IterSize iter_src) {
    auto iter_dst = order_.begin() + position_ * len_order_;
    if (!use_index_group_) {
      std::copy_n(iter_src,
                  len_order_,
                  iter_dst);
    } else {
      for (auto i : index_group_) {
        iter_dst = std::copy_n(iter_src + i * n_particles_,
                               n_particles_,
                               iter_dst);
      }
    }
  }
};

}
