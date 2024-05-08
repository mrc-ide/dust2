#pragma once

namespace dust2 {

template <typename T>
class dust_cpu {
public:
  using real_type = typename T::real_type;
  using rng_state_type = typename T::rng_state_type;
  using shared_state = typename T::shared_state;
  using internal_state = typename T::internal_state;
  using data_type = typename T::data_type;

  using rng_int_type = typename rng_state_type::int_type;

  dust_cpu(std::vector<shared_state> shared,
           std::vector<internal_state> internal,
           real_type time,
           real_type dt,
           size_t n_particles, // per group
           const std::vector<rng_int_type>& seed,
           bool deterministic) :
    n_state_(T::size(shared[0])),
    n_particles_(n_particles),
    n_groups_(shared.size()),
    n_particles_total_(n_particles_ * n_groups_),
    state_(n_state_ * n_particles_total_),
    state_next_(n_state_ * n_particles_total_),
    shared_(shared),
    internal_(internal),
    time_(time),
    dt_(dt),
    rng_(n_particles_total_, seed, deterministic) {
    // TODO: above, filter rng states need adding here too, or
    // somewhere at least (we might move the filter elsewhere though,
    // in which case that particular bit of weirdness goes away).

    // We don't check that the size is the same across all states;
    // this should be done by the caller (similarly, we don't check
    // that shared and internal have the same size).
  }

  auto run_steps(size_t n_steps) {
    // Ignore errors for now.
    real_type * state_data = state_.data();
    real_type * state_next_data = state_next_.data();

    // Later we parallelise this and track errors carefully.  We
    // should be able to parallelise with 'pragma omp parallel for
    // collapse(2)' here, but for MPI we might investigate other
    // approaches.
    for (size_t i = 0; i < n_groups_; ++i) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        run_particle(time_, dt_, n_steps,
                     shared_[i], internal_[i],
                     state_data + offset,
                     rng_.state(k),
                     state_next_data + offset);
      }
    }
    if (n_steps % 2 == 1) {
      std::swap(state_, state_next_);
    }
    // Time management here is going to require some effort once we
    // support interesting dt so that we always land on times with no
    // non-integer bits, but for now we require that dt is 1 so this
    // is easy.  We need this to hold within run_particle too, so it's
    // possible that's where the calculation here will be done.
    time_ = time_ + n_steps * dt_;
  }

  void set_state_initial() {
    real_type * state_data = state_.data();
    for (size_t i = 0; i < n_groups_; ++i) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        T::initial(time_, dt_,
                   shared_[i], internal_[i],
                   rng_.state(k),
                   state_data + offset);
      }
    }
  }

  template <typename It>
  void set_state(It it) {
    std::copy_n(it, state_.size(), state_.begin());
  }

  auto state() const {
    return state_;
  }

  // Fairly useless getter/setter - we might be better exposing time
  // directly as a field.  However, for the MPI and GPU version this
  // will almost certainly do something.
  auto time() const {
    return time_;
  }

  void set_time(real_type time) {
    time_ = time;
  }

  auto rng_state() { // TODO: should be const, error in mcstate2
    return rng_.export_state();
  }

  void set_rng_state(const std::vector<rng_int_type> & rng_state) {
    rng_.import_state();
  }

  template<typename Fn>
  void update_shared(size_t i, Fn fn) {
    // TODO: check that size was not modified, error if so (quite a
    // bit later).
    fn(shared_[i]);
  }

  auto n_groups() const {
    return n_groups_;
  }

private:
  size_t n_state_;
  size_t n_particles_;
  size_t n_groups_;
  size_t n_particles_total_;
  std::vector<real_type> state_;
  std::vector<real_type> state_next_;
  std::vector<shared_state> shared_;
  std::vector<internal_state> internal_;
  real_type time_;
  real_type dt_;
  mcstate::random::prng<rng_state_type> rng_;

  static void run_particle(real_type time, real_type dt, size_t n_steps,
                           shared_state shared, internal_state internal,
                           real_type * state, rng_state_type& rng_state,
                           real_type * state_next) {
    for (size_t i = 0; i < n_steps; ++i) {
      T::update(time + i * dt, dt, state, shared, internal, rng_state,
                state_next);
      std::swap(state, state_next);
    }
  }
};

}
