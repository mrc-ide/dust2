#pragma once

#include <algorithm>
#include <vector>
#include <dust2/continuous/control.hpp>
#include <dust2/continuous/solver.hpp>
#include <dust2/zero.hpp>
#include <mcstate/random/random.hpp>

namespace dust2 {

template <typename T>
class dust_continuous {
public:
  using system_type = T;
  using real_type = typename T::real_type;
  using rng_state_type = typename T::rng_state_type;
  using shared_state = typename T::shared_state;
  using internal_state = typename T::internal_state;
  using data_type = typename T::data_type;
  using rng_int_type = typename rng_state_type::int_type;

  dust_continuous(std::vector<shared_state> shared,
                  std::vector<internal_state> internal,
                  real_type time,
                  const ode::control<real_type> control, // in place of dt
                  size_t n_particles, // per group
                  const std::vector<rng_int_type>& seed,
                  bool deterministic) :
    n_state_(T::size_state(shared[0])),
    n_particles_(n_particles),
    n_groups_(shared.size()),
    n_particles_total_(n_particles_ * n_groups_),

    control_(control),

    state_(n_state_ * n_particles_total_),
    ode_internals_(n_particles_total_, n_state_),

    // For reordering to work:
    state_other_(n_state_ * n_particles_total_),
    ode_internals_other_(n_particles_total_, n_state_),

    shared_(shared),
    internal_(internal),

    time_(time),
    zero_every_(zero_every_vec<T>(shared_)),
    rng_(n_particles_total_, seed, deterministic),
    solver_(n_state_, control_) {
    // TODO: above, filter rng states need adding here too, or
    // somewhere at least (we might move the filter elsewhere though,
    // in which case that particular bit of weirdness goes away).

    // We don't check that the size is the same across all states;
    // this should be done by the caller (similarly, we don't check
    // that shared and internal have the same size).
  }

  auto run_to_time(real_type time) {
    real_type * state_data = state_.data();
    for (size_t i = 0; i < n_groups_; ++i) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        real_type * y = state_data + offset;
        solver_.run(time_, time, y, zero_every_[i],
                    ode_internals_[k],
                    rhs_(shared_[i], internal_[i]));
      }
    }
    time_ = time;
  }

  void set_state_initial() {
    real_type * state_data = state_.data();
    for (size_t i = 0; i < n_groups_; ++i) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        real_type * y = state_data + offset;
        T::initial(time_, shared_[i], internal_[i],
                   rng_.state(k), y);
        solver_.initialise(time_, y, ode_internals_[k],
                           rhs_(shared_[i], internal_[i]));
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
        const auto k = n_particles_ * i + j;
        const auto offset_read =
          i * offset_read_group + j * offset_read_particle;
        const auto offset_write = k * n_state_;
        real_type* y = state_data + offset_write;
        std::copy_n(iter + offset_read, n_state_, y);
        solver_.initialise(time_, y, ode_internals_[k],
                           rhs_(shared_[i], internal_[i]));
      }
    }
  }

  // To do this, we need a second copy of all internal state (so
  // y and internals).  But this is then a bit weird because we have
  // one block of memory that is really big (the 'y' bits) and passed
  // around as pointers and this other fragmented block.  The main aim
  // at this point should be just to get things implemented and then
  // we can consider efficiency and design later, but at least at the
  // moment we understand where the state actually comes from.
  //
  // Moving this onto the GPU, we'd want to do this differently, where
  // we hold these all as new elements.
  template <typename Iter>
  void reorder(Iter iter) {
    state_other_.resize(state_.size());
    for (size_t i = 0; i < n_groups_; ++i) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k_to = n_particles_ * i + j;
        const auto k_from = n_particles_ * i + *(iter + k_to);
        std::copy_n(state_.begin() + k_from * n_state_,
                    n_state_,
                    state_other_.begin() + k_to * n_state_);
        ode_internals_other_[k_to] = ode_internals_[k_from];
      }
    }
    std::swap(state_, state_other_);
    std::swap(ode_internals_, ode_internals_other_);
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

  auto rng_state() const {
    return rng_.export_state();
  }

  void set_rng_state(const std::vector<rng_int_type>& rng_state) {
    rng_.import_state(rng_state);
  }

  template<typename Fn>
  void update_shared(size_t i, Fn fn) {
    // TODO: check that size was not modified, error if so (quite a
    // bit later).
    fn(shared_[i]);
  }

  template <typename IterData, typename IterOutput>
  void compare_data(IterData data, IterOutput output) {
    const real_type * state_data = state_.data();
    for (size_t i = 0; i < n_groups_; ++i, ++data) {
      for (size_t j = 0; j < n_particles_; ++j, ++output) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        *output = T::compare_data(time_, state_data + offset, *data,
                                  shared_[i], internal_[i], rng_.state(k));
      }
    }
  }

  // This is just used for debugging
  const auto& ode_internals() const {
    return ode_internals_;
  }

private:
  size_t n_state_;
  size_t n_particles_;
  size_t n_groups_;
  size_t n_particles_total_;
  ode::control<real_type> control_;
  // Some more will be needed here to get history to work.  With
  // that, we'll need to hold something that will let us accumulate
  // (without allocation) into a large array which will be provided
  // to the integrator.  That could be owned by us, or by the
  // particle filter?
  std::vector<real_type> state_;
  std::vector<ode::internals<real_type>> ode_internals_;
  std::vector<real_type> state_other_;
  std::vector<ode::internals<real_type>> ode_internals_other_;
  std::vector<shared_state> shared_;
  std::vector<internal_state> internal_;
  real_type time_;
  std::vector<zero_every_type<real_type>> zero_every_;
  mcstate::random::prng<rng_state_type> rng_;
  ode::solver<real_type> solver_;

  static auto rhs_(const shared_state& shared, internal_state& internal) {
    return [&](real_type t, const real_type* y, real_type* dydt) {
             T::rhs(t, y, shared, internal, dydt);
           };
  }
};

}
