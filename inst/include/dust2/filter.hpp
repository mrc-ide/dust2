#pragma once

#include <dust2/filter_details.hpp>
#include <dust2/history.hpp>
#include <mcstate/random/random.hpp>

namespace dust2 {

template <typename T>
class filter {
public:
  using real_type = typename T::real_type;
  using data_type = typename T::data_type;
  using system_type = typename T::system_type;
  using rng_state_type = typename T::rng_state_type;
  using rng_int_type = typename rng_state_type::int_type;

  T sys;

  filter(T sys_,
         real_type time_start,
         std::vector<real_type> time,
         std::vector<data_type> data,
         std::vector<size_t> history_index,
         const std::vector<rng_int_type>& seed) :
    sys(sys_),
    time_start_(time_start),
    time_(time),
    data_(data),
    n_state_(sys.n_state()),
    n_particles_(sys.n_particles()),
    n_groups_(sys.n_groups()),
    rng_(n_groups_, seed, false),
    ll_(n_groups_ * n_particles_, 0),
    ll_step_(n_groups_ * n_particles_, 0),
    history_index_(history_index),
    history_(history_index_.size() > 0 ? history_index_.size() : n_state_,
             n_particles_, n_groups_, time_.size()),
    history_is_current_(false) {
    // TODO: duplicated with the above, can be done generically though
    // it's not a lot of code.
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

    // Just store this here; later once we have state to save we can
    // probably use that vector instead.
    std::vector<size_t> index(n_particles_ * n_groups_);

    const bool use_index = history_index_.size() > 0;

    auto it_data = data_.begin();
    for (size_t i = 0; i < n_times; ++i, it_data += n_groups_) {
      sys.run_steps(step_[i]);
      sys.compare_data(it_data, ll_step_.begin());

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

      sys.reorder(index.begin());

      if (save_history) {
        if (use_index) {
          history_.add_with_index(time_[i], sys.state().begin(), index.begin(),
                                  history_index_.begin(), n_state_);
        } else {
          history_.add(time_[i], sys.state().begin(), index.begin());
        }
      }
      // save snapshots (perhaps)
    }

    history_is_current_ = save_history;
  }

  template <typename It>
  void last_log_likelihood(It it) {
    std::copy_n(ll_.begin(), n_groups_, it);
  }

  auto& last_history() const {
    return history_;
  }

  bool last_history_is_current() const {
    return history_is_current_;
  }

  auto rng_state() const {
    return rng_.export_state();
  }

  auto set_rng_state(const std::vector<rng_int_type>& rng_state) {
    return rng_.import_state(rng_state);
  }

private:
  real_type time_start_;
  std::vector<real_type> time_;
  std::vector<size_t> step_;
  std::vector<data_type> data_;
  size_t n_state_;
  size_t n_particles_;
  size_t n_groups_;
  mcstate::random::prng<rng_state_type> rng_;
  std::vector<real_type> ll_;
  std::vector<real_type> ll_step_;
  std::vector<size_t> history_index_;
  history<real_type> history_;
  bool history_is_current_;
};

}
