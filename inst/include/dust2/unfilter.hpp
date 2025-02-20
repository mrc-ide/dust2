#pragma once

#include <dust2/adjoint_data.hpp>
#include <dust2/filter_details.hpp>
#include <dust2/trajectories.hpp>
#include <monty/random/random.hpp>

namespace dust2 {

template <typename T>
class unfilter {
public:
  using real_type = typename T::real_type;
  using system_type = typename T::system_type;
  using data_type = typename system_type::data_type;

  // We need to provide direct access to the system, because the user
  // will want to set parameters in, and pull out state, etc.
  T sys;

  unfilter(T sys_,
           real_type time_start,
           std::vector<real_type> time,
           std::vector<data_type> data) :
    sys(sys_),
    time_start_(time_start),
    time_(time),
    data_(data),
    n_state_(sys.n_state()),
    n_particles_(sys.n_particles()),
    n_groups_(sys.n_groups()),
    n_groups_data_(data_.size() / time.size()),
    ll_(n_particles_ * n_groups_, 0),
    ll_step_(n_particles_ * n_groups_, 0),
    trajectories_(n_state_, n_particles_, n_groups_, time_.size()),
    snapshots_(n_state_, n_particles_, n_groups_, 0),
    adjoint_(n_state_, n_particles_, n_groups_),
    trajectories_are_current_(n_particles_ * n_groups_),
    snapshots_are_current_(n_particles_ * n_groups_),
    adjoint_is_current_(n_particles_ * n_groups_),
    gradient_is_current_(false) {
  }

  void run(bool set_initial,
           bool save_trajectories,
           const std::vector<real_type>& save_snapshots,
           const std::vector<size_t>& index_state,
           const std::vector<size_t>& index_group) {
    const bool adjoint = false;
    reset(set_initial, save_trajectories, save_snapshots, index_state, index_group, adjoint);
    const auto n_times = time_.size();

    auto it_data = data_.begin();
    auto it_snapshots = save_snapshots.begin();
    for (size_t i = 0; i < n_times; ++i, it_data += n_groups_) {
      sys.run_to_time(time_[i], index_group);
      sys.compare_data(it_data, n_groups_data_, index_group, ll_step_.begin());
      for (size_t j = 0; j < ll_.size(); ++j) {
        ll_[j] += ll_step_[j];
      }
      if (save_trajectories) {
        trajectories_.add(time_[i], sys.state().begin());
      }
      if (it_snapshots != save_snapshots.end() && *it_snapshots == time_[i]) {
        snapshots_.add(time_[i], sys.state().begin());
        ++it_snapshots;
      }
    }

    if (save_trajectories) {
      for (auto i : index_group) {
        trajectories_are_current_[i] = true;
      }
    }
    if (!save_snapshots.empty()) {
      for (auto i : index_group) {
        snapshots_are_current_[i] = true;
      }
    }

    if (!tools::is_trivial_index(index_group, n_groups_)) {
      last_index_group_ = index_group;
    }
  }

  // This part here we can _always_ do, even if the system does not
  // actually support adjoint methods.  It should give exactly the
  // same answers as the normal version, at the cost of more memory.
  void run_adjoint(bool set_initial,
                   bool save_trajectories,
                   const std::vector<real_type> save_snapshots,
                   const std::vector<size_t>& index_state,
                   const std::vector<size_t>& index_group) {
    const bool adjoint = true;
    reset(set_initial, save_trajectories, save_snapshots, index_state, index_group, adjoint);

    // Run the entire forward time simulation
    sys.run_to_time(time_.back(), index_group, adjoint_.state(0));

    // Then all the data comparison in one pass.  This bit can
    // theoretically be parallelised but it's unlikely to be
    // important in most models.
    auto it_snapshots = save_snapshots.begin();
    const auto n_times = time_.size();
    for (size_t i = 0; i < n_times; ++i) {
      const auto state_i = adjoint_.state(i + 1);
      const auto data_i = data_.begin() + n_groups_ * i;
      sys.compare_data(data_i, n_groups_data_, state_i, index_group, ll_step_.begin());
      for (size_t j = 0; j < ll_.size(); ++j) {
        ll_[j] += ll_step_[j];
      }
      if (save_trajectories) {
        trajectories_.add(time_[i], state_i);
      }
      if (it_snapshots != save_snapshots.end() && *it_snapshots == time_[i]) {
        trajectories_.add(time_[i], state_i);
        ++it_snapshots;
      }
    }

    for (auto i : index_group) {
      adjoint_is_current_[i] = true;
      trajectories_are_current_[i] = save_trajectories;
      snapshots_are_current_[i] = !save_snapshots.empty();
    }

    if (!tools::is_trivial_index(index_group, n_groups_)) {
      last_index_group_ = index_group;
    }
  }

  auto& last_log_likelihood() const {
    return ll_;
  }

  auto& last_trajectories() const {
    return trajectories_;
  }

  auto& last_trajectories_are_current() const {
    return trajectories_are_current_;
  }

  auto& last_snapshots() const {
    return snapshots_;
  }

  auto& last_snapshots_are_current() const {
    return snapshots_are_current_;
  }

  auto& last_index_group() const {
    return last_index_group_.empty() ? sys.all_groups() : last_index_group_;
  }

  auto& adjoint_is_current() const {
    return adjoint_is_current_;
  }

  template <typename Iter>
  void last_gradient(Iter iter) {
    const auto& index_group = last_index_group();
    if (!gradient_is_current_) {
      compute_gradient_(index_group);
    }
    adjoint_.gradient(iter, index_group);
  }

private:
  real_type time_start_;
  std::vector<real_type> time_;
  std::vector<data_type> data_;
  size_t n_state_;
  size_t n_particles_;
  size_t n_groups_;
  size_t n_groups_data_;
  std::vector<real_type> ll_;
  std::vector<real_type> ll_step_;
  trajectories<real_type> trajectories_;
  trajectories<real_type> snapshots_;
  adjoint_data<real_type> adjoint_;
  std::vector<bool> trajectories_are_current_;
  std::vector<bool> snapshots_are_current_;
  std::vector<size_t> last_index_group_;
  std::vector<bool> adjoint_is_current_;
  bool gradient_is_current_;

  void reset(bool set_initial,
             bool save_trajectories,
             const std::vector<real_type>& save_snapshots,
             const std::vector<size_t>& index_state,
             const std::vector<size_t>& index_group,
             bool adjoint) {
    std::fill(trajectories_are_current_.begin(),
              trajectories_are_current_.end(),
              false);
    std::fill(snapshots_are_current_.begin(),
              snapshots_are_current_.end(),
              false);
    std::fill(adjoint_is_current_.begin(), adjoint_is_current_.end(), false);
    gradient_is_current_ = false;
    if (save_trajectories) {
      trajectories_.set_index_and_reset(index_state, index_group);
    }
    if (!save_snapshots.empty()) {
      snapshots_.set_n_times_and_reset(save_snapshots.size(), index_group);
    }
    if (adjoint) {
      adjoint_.init_history(time_start_, time_, sys.dt());
    }
    std::fill(ll_.begin(), ll_.end(), 0);
    sys.set_time(time_start_);
    if (set_initial) {
      sys.set_state_initial(sys.all_groups());
    }
    last_index_group_.clear();
  }

  void compute_gradient_(const std::vector<size_t>& index_group) {
    const auto n_times = time_.size();
    const auto n_adjoint = sys.n_adjoint();
    adjoint_.init_adjoint(n_adjoint);
    auto adjoint_curr = adjoint_.curr();
    auto adjoint_next = adjoint_.next();

    auto time = sys.time();

    for (size_t irev = 0; irev < n_times; ++irev) {
      const auto i = n_times - irev - 1;
      const auto time_i = i == 0 ? time_start_ : time_[i - 1];
      const auto state_i = adjoint_.state(i + 1);
      const auto data_i = data_.begin() + i * n_groups_;
      // Compare data
      sys.adjoint_compare_data(time, data_i, state_i,
                               n_adjoint, index_group,
                               adjoint_curr, adjoint_next);
      std::swap(adjoint_curr, adjoint_next);
      // Then run the system backwards from time => time_i
      const auto swap_adjoint = sys.adjoint_run_to_time(time, time_i, state_i,
                                                        n_adjoint, index_group,
                                                        adjoint_curr, adjoint_next);
      // Bookkeeping chore
      if (swap_adjoint) {
        std::swap(adjoint_curr, adjoint_next);
      }
      time = time_i;
    }

    // Initial conditions go right at the end, and are surprisingly
    // hard to work out.
    sys.adjoint_initial(time, adjoint_.state(0), n_adjoint, index_group,
                        adjoint_curr, adjoint_next);

    // At the end of the calculation, copy the final states so that
    // both copies within adjoint_ are the same - this means that the
    // gradient calculation will be correct.
    std::copy_n(adjoint_next, n_adjoint * n_groups_, adjoint_curr);

    gradient_is_current_ = true;
  }
};

}
