#pragma once

#include <vector>

namespace dust2 {

template <typename T>
class dust_cpu {
public:
  using model_type = T;
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

  auto run_to_time(real_type time) {
    const size_t n_steps = std::round(std::max(0.0, time - time_) / dt_);
    return run_steps(n_steps);
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
    time_ = time_ + n_steps * dt_;
  }

  // This is a special case for helping with the the adjoint model; we
  // run over all states and use the given space for storing as we go,
  // not the space in the model.  We *do* write out the final state
  // and time correctly, so that at the end of the simulation
  // everything is correct.
  auto run_steps(size_t n_steps, real_type * state_history) {
    const auto stride = n_state_ * n_particles_ * n_groups_;
    std::copy_n(state_.begin(), stride, state_history);
    for (size_t i = 0; i < n_groups_; ++i) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        auto state_model_ij = state_.data() + offset;
        auto state_history_ij = state_history + offset;
        run_particle(time_, dt_, n_steps, stride,
                     shared_[i], internal_[i],
                     state_model_ij, state_history_ij,
                     rng_.state(k));
      }
    }
    std::copy_n(state_history + stride * n_steps, stride, state_.begin());
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

  template <typename Iter>
  void set_state(Iter iter, bool recycle_particle, bool recycle_group) {
    const auto offset_read_group = recycle_group ? 0 :
      (n_state_ * (recycle_particle ? 1 : n_particles_));
    const auto offset_read_particle = recycle_particle ? 0 : n_state_;

    real_type * state_data = state_.data();
    for (size_t i = 0; i < n_groups_; ++i) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto offset_read =
          i * offset_read_group + j * offset_read_particle;
        const auto offset_write = (n_particles_ * i + j) * n_state_;
        std::copy_n(iter + offset_read,
                    n_state_,
                    state_data + offset_write);
      }
    }
  }

  template <typename Iter>
  void reorder(Iter iter) {
    for (size_t i = 0; i < n_groups_; ++i) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k_to = n_particles_ * i + j;
        const auto k_from = n_particles_ * i + *(iter + k_to);
        std::copy_n(state_.begin() + k_from * n_state_,
                    n_state_,
                    state_next_.begin() + k_to * n_state_);
      }
    }
    std::swap(state_, state_next_);
  }

  auto& state() const {
    return state_;
  }

  // Fairly useless getter/setter - we might be better exposing time
  // directly as a field.  However, for the MPI and GPU version this
  // will almost certainly do something.
  auto time() const {
    return time_;
  }

  auto dt() const {
    return dt_;
  }

  auto n_state() const {
    return n_state_;
  }

  auto n_particles() const {
    return n_particles_;
  }

  auto n_groups() const {
    return n_groups_;
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

  template <typename IterData, typename IterOutput>
  void compare_data(IterData data, IterOutput output) {
    compare_data(data, output, state_.data());
  }

  template <typename IterData, typename IterOutput>
  void compare_data(IterData data, IterOutput output, real_type * state) {
    for (size_t i = 0; i < n_groups_; ++i, ++data) {
      for (size_t j = 0; j < n_particles_; ++j, ++output) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        *output = T::compare_data(time_, dt_, state + offset, *data,
                                  shared_[i], internal_[i], rng_.state(k));
      }
    }
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
                           const shared_state& shared, internal_state& internal,
                           real_type * state, rng_state_type& rng_state,
                           real_type * state_next) {
    for (size_t i = 0; i < n_steps; ++i) {
      T::update(time + i * dt, dt, state, shared, internal, rng_state,
                state_next);
      std::swap(state, state_next);
    }
  }

  static void run_particle(real_type time, real_type dt,
                           size_t n_steps, size_t stride,
                           const shared_state& shared, internal_state& internal,
                           const real_type* state_model, real_type* state_history,
                           rng_state_type& rng_state) {
    auto state_curr = state_model;
    auto state_next = state_history + stride;
    for (size_t i = 0; i < n_steps; ++i, state_next += stride) {
      T::update(time + i * dt, dt, state_curr, shared, internal, rng_state,
                state_next);
      state_curr = state_next;
    }
  }
};

}
