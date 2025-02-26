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
    n_times_total_(n_times),
    len_state_(n_state_total_ * n_particles_ * n_groups_total_),
    len_order_(n_particles_ * n_groups_total_),
    position_state_(0),
    position_order_(0) {
  }

  // We might not actually store the index of time here, but something
  // else instead, because we don't really
  void set_index_and_reset(const std::vector<size_t>& index_state,
                           const std::vector<size_t>& index_group,
                           size_t n_times) {
    n_times_ = n_time;
    index_state_ = index_state;
    index_group_ = index_group;
    use_index_state_ = !tools::is_trivial_index(index_state_, n_state_total_);
    use_index_group_ = !tools::is_trivial_index(index_group_, n_groups_total_);
    n_state_ = use_index_state_ ? index_state.size() : n_state_total_;
    n_groups_ = use_index_group_ ? index_group.size() : n_groups_total_;
    len_state_ = n_state_ * n_particles_ * n_groups_;

    times_state_.resize(n_times_);
    state_.resize(len_state_ * n_times_);

    times_order_.resize(n_times_);
    order_.resize(len_order_ * n_times_total_);
    reorder_.resize(n_times_total_);

    position_state_ = 0;
    position_order_ = 0;
  }

  // in this case we are not reordering at all.
  template <typename IterReal>
  void add(real_type time, IterReal iter_state) {
    copy_state_(iter_state);
    times_state_[position_state_] = time;
    position_state_++;

    reorder_[position_order_] = false;
    position_order_++;
  }

  template <typename IterReal, typename IterSize>
  void add(real_type time, IterReal iter_state, IterSize iter_order) {
    copy_state_(iter_state);
    times_state_[position_state_] = time;
    position_state_++;

    // copy_order_(iter_order);
    // times_order_[position_order_] = time;
    // reorder_[position_order_] = true;
    // position_order_++;
    update_order(iter_order);
  }

  template <typename IterSize>
  void update_order(IterSize iter_order) {
    copy_order_(iter_order);
    times_order_[position_order_] = time;
    reorder_[position_order_] = true;
    position_order_++;
  }

  // These allow a consumer to allocate the right size structures for
  // time and state for the total that we've actually used.
  // auto size_time() const {
  //   return position_state_;
  // }

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
    return position_state_;
  }

  auto& index_group() const {
    return index_group_;
  }

  template <typename Iter>
  void export_time(Iter iter) const {
    std::copy_n(times_state_.begin(), position_state_, iter);
  }

  // Note that this uses 'select_particle' to avoid conflict with our
  // internal 'index_particle'. Also we are not extracting the
  // rectangle defined by the combination of the group/particle values
  // but pairs of (group,particle)'s.
  template <typename Iter>
  void export_state(Iter iter,
                    bool reorder,
                    const std::vector<size_t>& select_particle) const {
    reorder = reorder && n_particles_ > 1 && position_order_ > 0 &&
      tools::any(reorder_);

    auto use_select_particle = !select_particle.empty();

    if (!reorder && !use_select_particle) {
      // Optimised for simplest case of just dumping out everything
      std::copy_n(state_.begin(), position_state_ * len_state_, iter);
    } else if (!reorder) {
      for (size_t i = 0; i < position_state_; ++i) {
        const auto iter_state = state_.begin() + i * len_state_;
        for (size_t j = 0; j < n_groups_; ++j) {
          const auto offset_state = j * n_state_ * n_particles_;
          const auto offset_j = select_particle[j] * n_state_;
          iter = std::copy_n(iter_state + offset_state + offset_j,
                             n_state_,
                             iter);
        }
      }
    } else if (use_select_particle) {
    } else {
      for (size_t i = 0; i < n_groups_; ++j) {
        reorder_group_(i);
      }



    }
  }

private:
  size_t n_state_;
  size_t n_state_total_;
  size_t n_particles_;
  size_t n_groups_;
  size_t n_groups_total_;
  size_t n_times_;
  size_t n_times_total_;
  size_t len_state_; // length of an update to state
  size_t len_order_; // length of an update to order
  size_t position_state_;
  size_t position_order_;
  std::vector<real_type> times_state_;
  std::vector<real_type> times_order_;
  std::vector<real_type> state_;
  std::vector<size_t> order_;
  std::vector<bool> reorder_;
  std::vector<size_t> index_state_;
  std::vector<size_t> index_group_;
  std::vector<size_t> order_single_;
  std::vector<size_t> order_time_;
  bool std::vector<bool> reordered_;
  bool use_index_state_;
  bool use_index_group_;

  void apply_reorder() {
    if (!reordered_) {
      // Set up the empty order:
      for (size_t j = 0; j < n_groups; ++j) {
        for (size_t k = 0; k < n_particles; ++k) {
          order_time_[j * n_particles + k] = k;
        }
      }

      // We can do better here and just reorder state once we hit it I
      // think.
      for (size_t irev = 0; irev < n_times_; ++i) {
        const auto i = position_order_ - irev - 1;
        for (size_t j = 0; j < n_groups; ++j) {
          for (size_t k = 0; k < n_particles; ++k) {

          }
        }
      }
    }
  }


  template <typename Iter>
  void reorder_group_(typename std::vector<size_t>::iterator iter_order,
                      typename std::vector<size_t>::iterator index_particle) {
    for (size_t i = 0; i < n_particles_; ++i) {
      const auto index_particle_i = index_particle + i;
      *index_particle_i = *(iter_order + *index_particle_i);
    }

  }


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
