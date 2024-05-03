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

  dust_cpu(shared_state shared,
           internal_state internal,
           real_type time,
           real_type dt,
           size_t n_particles,
           const std::vector<rng_int_type>& seed,
           bool deterministic) :
    n_particles_(n_particles),
    n_state_(T::size(shared)),
    state_(n_state_ * n_particles_),
    state_next_(n_state_ * n_particles_),
    shared_(shared),
    internal_(internal),
    time_(time),
    dt_(dt),
    rng_(n_particles_, seed, deterministic) {
    // TODO: above, filter states need adding here too.
    if (dt != 1) {
      throw std::runtime_error("Requiring dt = 1 for now");
    }
  }

  auto run_steps(size_t n_steps) {
    // Ignore errors for now.
    real_type * state_data = state_.data();
    real_type * state_next_data = state_next_.data();
    // Later we parallelise this and track errors carefully.
    for (size_t i = 0; i < n_particles_; ++i) {
      const auto offset = i * n_state_;
      run_particle(time_, dt_, n_steps,
                   shared_, internal_,
                   state_data + offset,
                   rng_.state(i),
                   state_next_data + offset);
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
    for (size_t i = 0; i < n_particles_; ++i) {
      const auto offset = i * n_state_;
      T::initial(time_, dt_,
                 shared_, internal_,
                 rng_.state(i),
                 state_data + offset);
    }
  }

  template <typename It>
  void set_state(It it) {
    std::copy_n(it, state_.size(), state_.begin());
  }

  auto state() const {
    return state_;
  }

  auto time() const {
    return time_;
  }

  void set_time(real_type time) {
    // TODO: some validation needs to be done here to deal with
    // offsets relative to dt.
    time_ = time;
  }

  auto rng_state() { // TOOD: should be const, error in mcstate2
    return rng_.export_state();
  }

  void set_rng_state(const std::vector<rng_int_type> & rng_state) {
    rng_.import_state();
  }

private:
  size_t n_particles_;
  size_t n_state_;
  std::vector<real_type> state_;
  std::vector<real_type> state_next_;
  shared_state shared_;
  internal_state internal_;
  real_type time_;
  real_type dt_;
  mcstate::random::prng<rng_state_type> rng_;

  static void run_particle(real_type time, real_type dt, size_t n_steps,
                           shared_state shared, internal_state internal,
                           real_type * state, rng_state_type& rng_state,
                           real_type * state_next) {
    for (size_t i = 0; i < n_steps; ++i) {
      T::update(time, dt, state, shared, internal, rng_state, state_next);
      std::swap(state, state_next);
    }
  }
};


}
