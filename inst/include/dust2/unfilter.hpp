#pragma once

#include <dust2/adjoint_data.hpp>
#include <dust2/filter_details.hpp>
#include <dust2/history.hpp>
#include <mcstate/random/random.hpp>

namespace dust2 {

template <typename T>
class unfilter {
public:
  using real_type = typename T::real_type;
  using data_type = typename T::data_type;
  using system_type = typename T::system_type;

  // We need to provide direct access to the system, because the user
  // will want to set parameters in, and pull out state, etc.
  T sys;

  unfilter(T sys_,
           real_type time_start,
           std::vector<real_type> time,
           std::vector<data_type> data,
           std::vector<size_t> history_index_state) :
    sys(sys_),
    time_start_(time_start),
    time_(time),
    data_(data),
    n_state_(sys.n_state()),
    n_particles_(sys.n_particles()),
    n_groups_(sys.n_groups()),
    ll_(n_particles_ * n_groups_, 0),
    ll_step_(n_particles_ * n_groups_, 0),
    history_index_state_(history_index_state),
    history_(history_index_state_.size() > 0 ? history_index_state_.size() : n_state_,
             n_particles_, n_groups_, time_.size()),
    adjoint_(n_state_, n_particles_ * n_groups_),
    history_is_current_(n_particles_ * n_groups_),
    adjoint_is_current_(n_particles_ * n_groups_),
    gradient_is_current_(false) {
  }

  void run(bool set_initial, bool save_history) {
    run(set_initial, save_history, sys.all_groups());
  }

  void run(bool set_initial, bool save_history,
           const std::vector<size_t>& index_group) {
    reset(set_initial, save_history, /* adjoint = */ false);
    const auto n_times = time_.size();

    const bool use_index = history_index_state_.size() > 0;

    auto it_data = data_.begin();
    for (size_t i = 0; i < n_times; ++i, it_data += n_groups_) {
      sys.run_to_time(time_[i], index_group);
      sys.compare_data(it_data, ll_step_.begin(), index_group);
      for (size_t j = 0; j < ll_.size(); ++j) {
        ll_[j] += ll_step_[j];
      }
      if (save_history) {
        if (use_index) {
          history_.add_with_index(time_[i], sys.state().begin(),
                                  history_index_state_.begin(), n_state_,
                                  index_group);
        } else {
          history_.add(time_[i], sys.state().begin(),
                       index_group);
        }
      }
    }

    if (save_history) {
      for (auto i : index_group) {
        history_is_current_[i] = true;
      }
    }
  }

  // This part here we can _always_ do, even if the system does not
  // actually support adjoint methods.  It should give exactly the
  // same answers as the normal version, at the cost of more memory.
  void run_adjoint(bool set_initial, bool save_history) {
    run_adjoint(set_initial, save_history, sys.all_groups());
  }

  void run_adjoint(bool set_initial, bool save_history,
                   const std::vector<size_t>& index_group) {
    reset(set_initial, save_history, /* adjoint = */ true);

    // Run the entire forward time simulation
    auto state = adjoint_.state();
    sys.run_to_time(time_.back(), state, index_group);

    // Then all the data comparison in one pass.  This bit can
    // theoretically be parallelised but it's unlikely to be
    // important in most models.
    const bool use_index = history_index_state_.size() > 0;
    const auto dt = sys.dt();
    const auto n_times = time_.size();
    const auto stride_state = n_particles_ * n_groups_ * n_state_;
    for (size_t i = 0; i < n_times; ++i) {
      const size_t n_steps =
        std::round(std::max(0.0, time_[i] - time_start_) / dt);
      const auto state_i = state + n_steps * stride_state;
      const auto data_i = data_.begin() + n_groups_ * i;
      sys.compare_data(data_i, state_i, ll_step_.begin(), index_group);
      for (size_t j = 0; j < ll_.size(); ++j) {
        ll_[j] += ll_step_[j];
      }
      if (save_history) {
        if (use_index) {
          history_.add_with_index(time_[i], state_i,
                                  history_index_state_.begin(), n_state_,
                                  index_group);
        } else {
          history_.add(time_[i], state_i, index_group);
        }
      }
    }

    for (auto i : index_group) {
      adjoint_is_current_[i] = true;
      history_is_current_[i] = save_history;
    }
  }

  auto& last_log_likelihood() const {
    return ll_;
  }

  auto& last_history() const {
    return history_;
  }

  auto& last_history_is_current() const {
    return history_is_current_;
  }

  auto& adjoint_is_current() const {
    return adjoint_is_current_;
  }

  template <typename Iter>
  void last_gradient(Iter iter, const std::vector<size_t>& index_group) {
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
  std::vector<real_type> ll_;
  std::vector<real_type> ll_step_;
  std::vector<size_t> history_index_state_;
  history<real_type> history_;
  adjoint_data<real_type> adjoint_;
  std::vector<bool> history_is_current_;
  std::vector<bool> adjoint_is_current_;
  bool gradient_is_current_;

  void reset(bool set_initial, bool save_history, bool adjoint) {
    std::fill(history_is_current_.begin(), history_is_current_.end(), false);
    std::fill(adjoint_is_current_.begin(), adjoint_is_current_.end(), false);
    gradient_is_current_ = false;
    if (save_history) {
      history_.reset();
    }
    if (adjoint) {
      const auto dt = sys.dt();
      const size_t n_steps =
        std::round(std::max(0.0, time_.back() - time_start_) / dt);
      adjoint_.init_history(n_steps);
    }
    std::fill(ll_.begin(), ll_.end(), 0);
    sys.set_time(time_start_);
    if (set_initial) {
      sys.set_state_initial();
    }
  }

  void compute_gradient_() {
    compute_gradient_(sys.all_groups());
  }

  void compute_gradient_(const std::vector<size_t>& index_group) {
    const auto n_times = time_.size();
    const auto n_adjoint = sys.n_adjoint();
    adjoint_.init_adjoint(n_adjoint);
    auto adjoint_curr = adjoint_.curr();
    auto adjoint_next = adjoint_.next();

    const auto stride_state = n_particles_ * n_groups_ * n_state_;
    const auto state = adjoint_.state();
    const auto dt = sys.dt();
    auto time = sys.time();
    size_t position = std::round(std::max(0.0, time - time_start_) / dt);

    for (size_t irev = 0; irev < n_times; ++irev) {
      const auto i = n_times - irev - 1;
      const auto time_i = i == 0 ? time_start_ : time_[i - 1];
      const auto state_i = state + position * stride_state;
      const auto data_i = data_.begin() + i * n_groups_;
      // Compare data
      sys.adjoint_compare_data(time, data_i, state_i,
                               n_adjoint, adjoint_curr, adjoint_next,
                               index_group);
      std::swap(adjoint_curr, adjoint_next);
      // Then run the system backwards from time => time_i
      const auto n_steps = sys.adjoint_run_to_time(time, time_i, state_i,
                                                   n_adjoint,
                                                   adjoint_curr, adjoint_next,
                                                   index_group);
      // Bookkeeping chore
      if (n_steps % 2 == 1) {
        std::swap(adjoint_curr, adjoint_next);
      }
      position -= n_steps;
      time = time_i;
    }

    // Initial conditions go right at the end, and are surprisingly
    // hard to work out.
    sys.adjoint_initial(time, state, n_adjoint, adjoint_curr, adjoint_next,
                        index_group);

    // At the end of the calculation, copy the final states so that
    // both copies within adjoint_ are the same - this means that the
    // gradient calculation will be correct.
    std::copy_n(adjoint_next, n_adjoint * n_groups_, adjoint_curr);

    gradient_is_current_ = true;
  }
};

}
