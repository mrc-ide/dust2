#pragma once

#include <dust2/cpu.hpp>

namespace dust2 {

template <typename T>
class unfilter {
public:
  using real_type = typename T::real_type;
  using data_type = typename T::data_type;

  // We need to provide direct access to the model, because the user
  // will want to set parameters in, and pull out state, etc.
  dust_cpu<T> model;

  unfilter(dust_cpu<T> model_,
           real_type time_start,
           std::vector<real_type> time,
           std::vector<data_type> data) :
    model(model_),
    time_start_(time_start),
    time_(time),
    data_(data),
    n_particles_(model.n_particles()),
    n_groups_(model.n_groups()),
    ll_(n_particles_ * n_groups_, 0),
    ll_step_(n_particles_ * n_groups_, 0) {
    const auto dt = model_.dt();
    for (size_t i = 0; i < time_.size(); i++) {
      const auto t0 = i == 0 ? time_start_ : time_[i - 1];
      const auto t1 = time_[i];
      step_.push_back(static_cast<size_t>(std::round((t1 - t0) / dt)));
    }
  }

  void run(bool set_initial) {
    const auto n_times = step_.size();

    model.set_time(time_start_);
    if (set_initial) {
      model.set_state_initial();
    }
    std::fill(ll_.begin(), ll_.end(), 0);

    auto it_data = data_.begin();
    for (size_t i = 0; i < n_times; ++i, it_data += n_groups_) {
      model.run_steps(step_[i]); // just compute this at point of use?
      model.compare_data(it_data, ll_step_.begin());
      for (size_t j = 0; j < ll_.size(); ++j) {
        ll_[j] += ll_step_[j];
      }
    }
  }

  template <typename Iter>
  void last_log_likelihood(Iter iter) {
    std::copy(ll_.begin(), ll_.end(), iter);
  }

private:
  real_type time_start_;
  std::vector<real_type> time_;
  std::vector<size_t> step_;
  std::vector<data_type> data_;
  size_t n_particles_;
  size_t n_groups_;
  std::vector<real_type> ll_;
  std::vector<real_type> ll_step_;
};

}
