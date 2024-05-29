#pragma once

#include <dust2/cpu.hpp>
#include <dust2/filter_details.hpp>
#include <mcstate/random/random.hpp>

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

template <typename T>
class filter {
public:
  using real_type = typename T::real_type;
  using data_type = typename T::data_type;
  using rng_state_type = typename T::rng_state_type;
  using rng_int_type = typename rng_state_type::int_type;

  dust_cpu<T> model;

  filter(dust_cpu<T> model_,
         real_type time_start,
         std::vector<real_type> time,
         std::vector<data_type> data,
         const std::vector<rng_int_type>& seed) :
    model(model_),
    time_start_(time_start),
    time_(time),
    data_(data),
    n_particles_(model.n_particles()),
    n_groups_(model.n_groups()),
    rng_(n_groups_, seed, false),
    ll_(n_groups_ * n_particles_, 0),
    ll_step_(n_groups_ * n_particles_, 0) {
    // TODO: duplicated with the above, can be done generically though
    // it's not a lot of code.
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

    // Just store this here; later once we have state to save we can
    // probably use that vector instead.
    std::vector<size_t> index(n_particles_ * n_groups_);

    auto it_data = data_.begin();
    for (size_t i = 0; i < n_times; ++i, it_data += n_groups_) {
      model.run_steps(step_[i]);
      model.compare_data(it_data, ll_step_.begin());

      for (size_t i = 0; i < n_groups_; ++i) {
        const auto offset = i * n_particles_;
        const auto w = ll_step_.begin() + offset;
        ll_[i] += details::scale_log_weights<real_type>(n_particles_, w);
      }

      // early exit here, once enabled, setting and checking a
      // threshhold log likelihood below which we're uninterested;
      // this requires some fiddling.

      // This can be parallelised across groups
      for (size_t i = 0; i < n_groups_; ++i) {
        const auto offset = i * n_particles_;
        const auto w = ll_step_.begin() + offset;
        const auto idx = index.begin() + offset;
        const auto u = mcstate::random::random_real<real_type>(rng_.state(i));
        details::resample_weight(n_particles_, w, u, idx);
      }

      model.reorder(index.begin());

      // save trajectories (perhaps)
      // save snapshots (perhaps)
    }
  }

  template <typename It>
  void last_log_likelihood(It it) {
    std::copy_n(ll_.begin(), n_groups_, it);
  }

  auto rng_state() { // TODO: should be const, error in mcstate2
    return rng_.export_state();
  }

private:
  real_type time_start_;
  std::vector<real_type> time_;
  std::vector<size_t> step_;
  std::vector<data_type> data_;
  size_t n_particles_;
  size_t n_groups_;
  mcstate::random::prng<rng_state_type> rng_;
  std::vector<real_type> ll_;
  std::vector<real_type> ll_step_;
};

}
