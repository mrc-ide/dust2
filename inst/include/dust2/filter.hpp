#pragma once

#include <dust2/filter_details.hpp>
#include <dust2/history.hpp>
#include <monty/random/random.hpp>

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
         std::vector<size_t> history_index_state,
         const std::vector<rng_int_type>& seed) :
    sys(sys_),
    time_start_(time_start),
    time_(time),
    data_(data),
    n_state_(sys.n_state()),
    n_particles_(sys.n_particles()),
    n_groups_(sys.n_groups()),
    rng_(n_groups_, seed, false),
    ll_(n_groups_, 0),
    ll_step_(n_groups_ * n_particles_, 0),
    history_index_state_(history_index_state),
    history_(history_index_state_.size() > 0 ? history_index_state_.size() : n_state_,
             n_particles_, n_groups_, time_.size()),
    history_is_current_(n_groups_, false) {
  }

  void run(bool set_initial, bool save_history,
           const std::vector<size_t>& index_group) {
    reset(set_initial, save_history);

    // Just store this here; later once we have state to save we can
    // probably use that vector instead.
    std::vector<size_t> index(n_particles_ * n_groups_);

    const bool use_index = history_index_state_.size() > 0;

    auto it_data = data_.begin();
    for (auto time : time_) {
      sys.run_to_time(time, index_group);
      sys.compare_data(it_data, index_group, ll_step_.begin());
      it_data += n_groups_;

      // This looks like it's just not capturing the jump properly
      // perhaps?  Better would be to return nan from the comparison
      // and cope with that accordingly - all nan is missing, some nan
      // is failure?
      for (auto i : index_group) {
        const auto offset = i * n_particles_;
        const auto w = ll_step_.begin() + offset;
	const auto idx = index.begin() + offset;

	// Only reorder a group if any value has done a comparison (in
	// the case where all ll values are exactly zero we take this
	// as the case where data was all empty).
	const bool reorder_this_group =
	  std::any_of(w, w + n_particles_, [](auto v) { return v != 0; });
	if (reorder_this_group) {
	  ll_[i] += details::scale_log_weights<real_type>(n_particles_, w);
	  const auto u = monty::random::random_real<real_type>(rng_.state(i));
	  details::resample_weight(n_particles_, w, u, idx);
	} else {
	  for (size_t j = 0; j < n_particles_; ++j) {
	    *(idx + j) = j;
	  }
	}
      }

      sys.reorder(index.begin(), index_group);

      if (save_history) {
        if (use_index) {
          history_.add_with_index(time, sys.state().begin(), index.begin(),
                                  history_index_state_.begin(), n_state_,
                                  index_group);
        } else {
          history_.add(time, sys.state().begin(), index.begin(),
                       index_group);
        }
      }
      // early exit (perhaps)
      // save snapshots (perhaps)
    }

    if (save_history) {
      for (auto i : index_group) {
        history_is_current_[i] = true;
      }
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

  auto rng_state() const {
    return rng_.export_state();
  }

  auto set_rng_state(const std::vector<rng_int_type>& rng_state) {
    return rng_.import_state(rng_state);
  }

  auto select_random_particle(size_t group) {
    const auto u = mcstate::random::random_real<real_type>(rng_.state(group));
    return static_cast<size_t>(std::floor(u * n_particles_));
  }

private:
  real_type time_start_;
  std::vector<real_type> time_;
  std::vector<data_type> data_;
  size_t n_state_;
  size_t n_particles_;
  size_t n_groups_;
  monty::random::prng<rng_state_type> rng_;
  std::vector<real_type> ll_;
  std::vector<real_type> ll_step_;
  std::vector<size_t> history_index_state_;
  history<real_type> history_;
  std::vector<bool> history_is_current_;

  void reset(bool set_initial, bool save_history) {
    std::fill(history_is_current_.begin(), history_is_current_.end(), false);
    if (save_history) {
      history_.reset();
    }
    std::fill(ll_.begin(), ll_.end(), 0);
    sys.set_time(time_start_);
    if (set_initial) {
      sys.set_state_initial(sys.all_groups());
    }
  }
};

}
