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
    n_snapshots_(0),
    len_state_(n_state_total_ * n_particles_ * n_groups_total_),
    len_order_(n_particles_ * n_groups_total_),
    position_state_(0),
    position_snapshot_(0),
    position_order_(0),
    save_state_(false) {
  }

  // We might not actually store the index of time here, but something
  // else instead, because we don't really
  void set_index_and_reset(const std::vector<size_t>& index_state,
                           const std::vector<size_t>& index_group,
                           bool save_state,
                           const std::vector<real_type>& times_snapshot) {
    save_state_ = save_state,
    index_state_ = index_state;
    index_group_ = index_group;
    times_snapshot_ = times_snapshot;
    use_index_state_ = !tools::is_trivial_index(index_state_, n_state_total_);
    use_index_group_ = !tools::is_trivial_index(index_group_, n_groups_total_);
    n_state_ = use_index_state_ ? index_state.size() : n_state_total_;
    n_groups_ = use_index_group_ ? index_group.size() : n_groups_total_;
    n_snapshots_ = times_snapshot_.size();
    len_state_ = n_state_ * n_particles_ * n_groups_;

    times_state_.resize(n_times_);
    times_order_.resize(n_times_);
    state_.resize(len_state_ * n_times_);
    snapshots_.resize(n_state_total_ * n_particles_ * n_groups_ * n_snapshots_);

    order_.resize(len_order_ * n_times_);
    reorder_.resize(n_times_);

    position_state_ = 0;
    position_order_ = 0;
    position_snapshot_ = 0;
  }

  // in this case we are not reordering at all.
  template <typename IterReal>
  void add(real_type time, IterReal iter_state) {
    add_state(time, iter_state);

    reorder_[position_order_] = false;
    position_order_++;
  }

  template <typename IterReal, typename IterSize>
  void add(real_type time, IterReal iter_state, IterSize iter_order) {
    add_state(time, iter_state);

    copy_order_(iter_order);
    times_order_[position_order_] = time;
    reorder_[position_order_] = true;
    position_order_++;
  }

  template <typename IterReal>
  void add_state(real_type time, IterReal iter_state) {
    if (save_state_) {
      copy_state_(iter_state);
      times_state_[position_state_] = time;
      position_state_++;
    }

    const bool save_snapshot =
      position_snapshot_ < n_snapshots_ &&
      time == times_snapshot_[position_snapshot_];
    if (save_snapshot) {
      copy_snapshot_(iter_state);
      position_snapshot_++;
    }
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
    return position_state_;
  }

  auto n_snapshots() const {
    return position_snapshot_;
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
    } else {
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
      for (size_t irev = 0; irev < position_state_; ++irev) {
        const auto i = position_state_ - irev - 1;
        const auto iter_order = order_.begin() + i * len_order_;
        const auto iter_state = state_.begin() + i * len_state_;
        for (size_t j = 0; j < n_groups_; ++j) {
          const auto offset_state = j * n_state_ * n_particles_;
          const auto offset_index = j * n_particles_;
          const auto offset_output =
            i * len_state_out + j * n_state_ * n_particles_out;
          if (use_select_particle) {
            reorder_state_single_(iter_state + offset_state,
                            iter_order + offset_index,
                            reorder_[i],
                            iter + offset_output,
                            index_particle[j]);
          } else {
            reorder_state_group_(iter_state + offset_state,
                           iter_order + offset_index,
                           reorder_[i],
                           iter + offset_output,
                           index_particle.begin() + offset_index);
          }
        }
      }
    }
  }

  // There is some opportunity to harmonise with above, but we'll get
  // things working before trying that.
  template <typename Iter>
  void export_snapshots(Iter iter,
                        bool reorder,
                        const std::vector<size_t>& select_particle) const {
    reorder = reorder && n_particles_ > 1 && position_order_ > 0 &&
      tools::any(reorder_);
    auto use_select_particle = !select_particle.empty();

    if (!reorder && !use_select_particle) {
      // Optimised for simplest case of just dumping out everything
      const auto len = n_state_total_ * n_particles_ * n_groups_total_;
      std::copy_n(snapshots_.begin(), position_snapshot_ * len, iter);
    } else if (!reorder) {
      const auto len = n_state_total_ * n_particles_ * n_groups_total_;
      for (size_t i = 0; i < position_snapshot_; ++i) {
        const auto iter_snapshot = snapshots_.begin() + i * len;
        for (size_t j = 0; j < n_groups_; ++j) {
          const auto offset_snapshot = j * n_state_total_ * n_particles_;
          const auto offset_j = select_particle[j] * n_state_total_;
          iter = std::copy_n(iter_snapshot + offset_snapshot + offset_j,
                             n_state_total_,
                             iter);
        }
      }
    } else {
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

      const auto len_state_out = n_state_total_ * n_particles_out * n_groups_;
      size_t i_snapshot = n_snapshots_ - 1;
      real_type time_next_snapshot = times_snapshot_[i_snapshot];

      for (size_t irev = 0; irev < position_order_; ++irev) {
        const auto i = position_order_ - irev - 1;
        const auto time = times_order_[i];

        if (time == time_next_snapshot) {
          const auto iter_src = snapshots_.begin() +
            i_snapshot * n_state_total_ * n_particles_ * n_groups_total_;
          const auto iter_dest = iter + i_snapshot * len_state_out;
          for (size_t j = 0; j < n_groups_; ++j) {
            if (use_select_particle) {
              const auto index_particle_j = index_particle[j];
              auto iter_src_j = iter_src +
                j * n_state_total_ * n_particles_ +
                index_particle_j * n_state_total_;
              auto iter_dest_j = iter_dest + j * n_state_total_;
              std::copy_n(iter_src_j, n_state_total_, iter_dest_j);
            } else {
              auto iter_dest_j = iter_dest + j * n_state_total_;
              for (size_t k = 0; k < n_particles_; ++k) {
                auto iter_src_k = iter_src +
                  j * n_state_total_ * n_particles_ +
                  index_particle[j * n_particles_ + k] * n_state_total_;
                auto iter_dest_k = iter_dest_j + k * n_state_total_;
                std::copy_n(iter_src_k, n_state_total_, iter_dest_k);
              }
            }
          }

          if (i_snapshot == 0) {
            break;
          }
          --i_snapshot;
        }

        if (reorder_[i]) {
          const auto iter_order = order_.begin() + i * len_order_;
          for (size_t j = 0; j < n_groups_; ++j) {
            const auto iter_order_j = iter_order + j * n_particles_;
            if (use_select_particle) {
              index_particle[j] = *(iter_order_j + index_particle[j]);
            } else {
              const auto index_particle_j = index_particle.begin() + j * n_particles_;
              for (size_t k = 0; k < n_groups_; ++k) {
                const auto index_particle_k = index_particle_j + k;
                *index_particle_k = *(iter_order_j + *index_particle_k);
              }
            }
          }
        }
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
  size_t n_snapshots_;
  size_t len_state_; // length of an update to state
  size_t len_order_; // length of an update to order
  size_t position_state_;
  size_t position_snapshot_;
  size_t position_order_;
  std::vector<real_type> times_state_;
  std::vector<real_type> times_snapshot_;
  std::vector<real_type> times_order_;
  std::vector<real_type> state_;
  std::vector<real_type> snapshots_;
  std::vector<size_t> order_;
  std::vector<bool> reorder_;
  std::vector<size_t> index_state_;
  std::vector<size_t> index_group_;
  std::vector<size_t> order_single_;
  std::vector<size_t> order_time_;
  std::vector<bool> reordered_;
  bool save_state_;
  bool use_index_state_;
  bool use_index_group_;

  // Reference implementation for this is mcstate::history_single and
  // mcstate::history_multiple
  template <typename Iter>
  void reorder_state_group_(typename std::vector<real_type>::const_iterator iter_state,
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
  void reorder_state_single_(typename std::vector<real_type>::const_iterator iter_state,
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
    auto iter_dst = state_.begin() + position_state_ * len_state_;
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

  template <typename IterReal>
  void copy_snapshot_(IterReal iter_src) {
    const auto len = n_state_total_ * n_particles_ * n_groups_;
    auto iter_dst = snapshots_.begin() + position_snapshot_ * len;
    if (!use_index_group_) {
      std::copy_n(iter_src, len, iter_dst);
    } else {
      const auto len_i = n_state_total_ * n_particles_;
      for (auto i : index_group_) {
        const auto iter_src_i = iter_src + i * len_i;
        iter_dst = std::copy_n(iter_src_i, len_i, iter_dst);
      }
    }
  }

  template <typename IterSize>
  void copy_order_(IterSize iter_src) {
    auto iter_dst = order_.begin() + position_order_ * len_order_;
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
