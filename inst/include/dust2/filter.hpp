#pragma once

#include <dust2/cpu.hpp>
#include <dust2/filter_details.hpp>
#include <dust2/history.hpp>
#include <dust2/adjoint.hpp>
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
           std::vector<data_type> data,
           std::vector<size_t> history_index) :
    model(model_),
    time_start_(time_start),
    time_(time),
    data_(data),
    n_state_(model.n_state()),
    n_particles_(model.n_particles()),
    n_groups_(model.n_groups()),
    ll_(n_particles_ * n_groups_, 0),
    ll_step_(n_particles_ * n_groups_, 0),
    history_index_(history_index),
    history_(history_index_.size() > 0 ? history_index_.size() : n_state_,
             n_particles_, n_groups_, time_.size()),
    adjoint_(n_state_, n_particles_ * n_groups_, time_.size()),
    history_is_current_(false),
    adjoint_is_current_(false),
    gradient_is_current_(false) {
    const auto dt = model_.dt();
    for (size_t i = 0; i < time_.size(); i++) {
      const auto t0 = i == 0 ? time_start_ : time_[i - 1];
      const auto t1 = time_[i];
      const size_t n = std::round((t1 - t0) / dt);
      step_.push_back(n);
      step_tot_.push_back(i == 0 ? n : step_tot_[i - 1] + n);
    }
  }

  void run(bool set_initial, bool save_history) {
    reset(save_history, false);
    const auto n_times = step_.size();

    model.set_time(time_start_);
    if (set_initial) {
      model.set_state_initial();
    }
    std::fill(ll_.begin(), ll_.end(), 0);

    const bool use_index = history_index_.size() > 0;

    auto it_data = data_.begin();
    for (size_t i = 0; i < n_times; ++i, it_data += n_groups_) {
      model.run_steps(step_[i]); // just compute this at point of use?
      model.compare_data(it_data, ll_step_.begin());
      for (size_t j = 0; j < ll_.size(); ++j) {
        ll_[j] += ll_step_[j];
      }
      if (save_history) {
        if (use_index) {
          history_.add_with_index(time_[i], model.state().begin(),
                                  history_index_.begin(), n_state_);
        } else {
          history_.add(time_[i], model.state().begin());
        }
      }
    }

    history_is_current_ = save_history;
  }

  // This part here we can _always_ do, even if the model does not
  // actually support adjoint methods.  It should give exactly the
  // same answers as the normal version, at the cost of more memory.
  void run_adjoint(bool set_initial, bool save_history) {
    reset(save_history, true);
    const auto n_times = step_.size();
    if (set_initial) {
      model.set_state_initial();
    }
    std::fill(ll_.begin(), ll_.end(), 0);

    // Run the entire forward time simulation
    auto state = adjoint_.state();
    model.run_steps(step_tot_.back(), state);

    // Consider dumping out a fraction of history here into our normal
    // history saving; that's just a copy really.
    if (save_history) {
      throw std::runtime_error("not yet implemented");
    }

    // Then all the data comparison in one pass.
    const auto stride_state = n_particles_ * n_groups_ * n_state_;
    auto it_data = data_.begin();
    for (size_t i = 0; i < n_times; ++i, it_data += n_groups_) {
      state += step_[i] * stride_state;
      model.compare_data(it_data, ll_step_.begin(), state);
      for (size_t j = 0; j < ll_.size(); ++j) {
        ll_[j] += ll_step_[j];
      }
    }

    adjoint_is_current_ = true;
    history_is_current_ = save_history;
  }

  template <typename Iter>
  void last_log_likelihood(Iter iter) {
    std::copy(ll_.begin(), ll_.end(), iter);
  }

  auto& last_history() const {
    // In the case where adjoint_is_current_ &&
    // !history_is_current_, we can fairly efficiently copy the
    // history over and then return, though that means that this is no
    // longer a const method (but the return value should still be
    // marked as such).  If we do that then the test below should be
    // ||'d with adjoint_is_current_.
    return history_;
  }

  bool last_history_is_current() const {
    return history_is_current_;
  }

  template <typename Iter>
  void last_gradient(Iter iter) {
    if (!gradient_is_current_) {
      compute_gradient_();
    }
    adjoint_.gradient(iter);
  }

private:
  real_type time_start_;
  std::vector<real_type> time_;
  std::vector<size_t> step_;
  std::vector<size_t> step_tot_;
  std::vector<data_type> data_;
  size_t n_state_;
  size_t n_particles_;
  size_t n_groups_;
  std::vector<real_type> ll_;
  std::vector<real_type> ll_step_;
  std::vector<size_t> history_index_;
  history<real_type> history_;
  adjoint_data<real_type> adjoint_;
  bool history_is_current_;
  bool adjoint_is_current_;
  bool gradient_is_current_;

  void reset(bool save_history, bool adjoint) {
    history_is_current_ = false;
    adjoint_is_current_ = false;
    gradient_is_current_ = false;
    if (save_history) {
      history_.reset();
    }
    if (adjoint) {
      adjoint_.init_history(step_tot_.back());
    }
  }

  void compute_gradient_() {
    // Probably we just need this on the R interface? We don't need to
    // expose the gradient_is_current part though.
    if (!adjoint_is_current_) {
      throw std::runtime_error("adjoint is not current");
    }

    const auto n_times = time_.size();
    const auto n_adjoint = T::size_adjoint();
    const auto adjoint_curr = adjoint_.curr();
    const auto adjoint_next = adjoint_.next();
    auto it_data = data_.last();
    
    const auto stride_state = n_particles_ * n_groups_ * n_state_;
    const auto state = adjoint_.state();
    const auto dt = model.dt();

    // We do need the time here, and there are a couple of ways of
    // getting it.
    for (size_t irev = 0; irev < n_times; ++irev, it_data -= n_groups_) {
      const auto i = n_times - irev - 1;
      const auto time = time_start_ + step_tot_[i] * dt;
      const auto n_steps = step_[i];
      const auto state_i = state + step_tot_[i] * stride_state;
      // Compare data
      model.adjoint_compare_data(it_data, state_i, adjoint_curr, adjoint_next);
      std::swap(adjoint_curr, adjoint_next);
      // Then run the model backwards
      model.adjoint_run_steps(n_steps, time,
                              state, adjoint_curr, adjoint_next);
      // Bookkeeping chore
      if (n_steps % 2 == 1) {
        std::swap(adjoint_curr, adjoint_next);
      }
    }

    // Initial conditions go right at the end, and are surprisingly
    // hard to work out.
    model.adjoint_initial(time_start_, state, adjoint_curr, adjoint_next);
    std::swap(adjoint_curr, adjoint_next);
    gradient_is_current_ = true;
  }
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
         std::vector<size_t> history_index,
         const std::vector<rng_int_type>& seed) :
    model(model_),
    time_start_(time_start),
    time_(time),
    data_(data),
    n_state_(model.n_state()),
    n_particles_(model.n_particles()),
    n_groups_(model.n_groups()),
    rng_(n_groups_, seed, false),
    ll_(n_groups_ * n_particles_, 0),
    ll_step_(n_groups_ * n_particles_, 0),
    history_index_(history_index),
    history_(history_index_.size() > 0 ? history_index_.size() : n_state_,
             n_particles_, n_groups_, time_.size()),
    history_is_current_(false) {
    // TODO: duplicated with the above, can be done generically though
    // it's not a lot of code.
    const auto dt = model_.dt();
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

    model.set_time(time_start_);
    if (set_initial) {
      model.set_state_initial();
    }
    std::fill(ll_.begin(), ll_.end(), 0);

    // Just store this here; later once we have state to save we can
    // probably use that vector instead.
    std::vector<size_t> index(n_particles_ * n_groups_);

    const bool use_index = history_index_.size() > 0;

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

      if (save_history) {
        if (use_index) {
          history_.add_with_index(time_[i], model.state().begin(), index.begin(),
                                  history_index_.begin(), n_state_);
        } else {
          history_.add(time_[i], model.state().begin(), index.begin());
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

  auto rng_state() { // TODO: should be const, error in mcstate2
    return rng_.export_state();
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
