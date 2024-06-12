#pragma once

#include <dust2/discrete/system.hpp>
#include <dust2/filter_details.hpp>
#include <dust2/history.hpp>
#include <mcstate/random/random.hpp>

namespace dust2 {

template <typename T>
class unfilter {
public:
  using real_type = typename T::real_type;
  using data_type = typename T::data_type;

  // We need to provide direct access to the system, because the user
  // will want to set parameters in, and pull out state, etc.
  dust_discrete<T> sys;

  unfilter(dust_discrete<T> sys_,
           real_type time_start,
           std::vector<real_type> time,
           std::vector<data_type> data,
           std::vector<size_t> history_index) :
    sys(sys_),
    time_start_(time_start),
    time_(time),
    data_(data),
    n_state_(sys.n_state()),
    n_particles_(sys.n_particles()),
    n_groups_(sys.n_groups()),
    ll_(n_particles_ * n_groups_, 0),
    ll_step_(n_particles_ * n_groups_, 0),
    history_index_(history_index),
    history_(history_index_.size() > 0 ? history_index_.size() : n_state_,
             n_particles_, n_groups_, time_.size()),
    history_is_current_(false) {
    const auto dt = sys_.dt();
    for (size_t i = 0; i < time_.size(); i++) {
      const auto t0 = i == 0 ? time_start_ : time_[i - 1];
      const auto t1 = time_[i];
      step_.push_back(static_cast<size_t>(std::round((t1 - t0) / dt)));
    }
  }

  void run(bool set_initial, bool save_history) {
    history_is_current_ = false;
    if (save_history) {
      history_.reset();
    }
    const auto n_times = step_.size();

    sys.set_time(time_start_);
    if (set_initial) {
      sys.set_state_initial();
    }
    std::fill(ll_.begin(), ll_.end(), 0);

    const bool use_index = history_index_.size() > 0;

    auto it_data = data_.begin();
    for (size_t i = 0; i < n_times; ++i, it_data += n_groups_) {
      sys.run_steps(step_[i]); // just compute this at point of use?
      sys.compare_data(it_data, ll_step_.begin());
      for (size_t j = 0; j < ll_.size(); ++j) {
        ll_[j] += ll_step_[j];
      }
      if (save_history) {
        if (use_index) {
          history_.add_with_index(time_[i], sys.state().begin(),
                                  history_index_.begin(), n_state_);
        } else {
          history_.add(time_[i], sys.state().begin());
        }
      }
    }

    history_is_current_ = save_history;
  }

  template <typename Iter>
  void last_log_likelihood(Iter iter) {
    std::copy(ll_.begin(), ll_.end(), iter);
  }


  auto& last_history() const {
    return history_;
  }

  bool last_history_is_current() const {
    return history_is_current_;
  }

private:
  real_type time_start_;
  std::vector<real_type> time_;
  std::vector<size_t> step_;
  std::vector<data_type> data_;
  size_t n_state_;
  size_t n_particles_;
  size_t n_groups_;
  std::vector<real_type> ll_;
  std::vector<real_type> ll_step_;
  std::vector<size_t> history_index_;
  history<real_type> history_;
  bool history_is_current_;
};

}
