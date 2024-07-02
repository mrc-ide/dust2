#pragma once

#include <map>
#include <vector>
#include <dust2/zero.hpp>
#include <mcstate/random/random.hpp>

namespace dust2 {

template <typename T>
class dust_discrete {
public:
  using system_type = T;
  using real_type = typename T::real_type;
  using rng_state_type = typename T::rng_state_type;
  using shared_state = typename T::shared_state;
  using internal_state = typename T::internal_state;
  using data_type = typename T::data_type;

  using rng_int_type = typename rng_state_type::int_type;

  dust_discrete(std::vector<shared_state> shared,
                std::vector<internal_state> internal,
                real_type time,
                real_type dt,
                size_t n_particles, // per group
                const std::vector<rng_int_type>& seed,
                bool deterministic) :
    n_state_(T::size_state(shared[0])),
    n_particles_(n_particles),
    n_groups_(shared.size()),
    n_particles_total_(n_particles_ * n_groups_),
    state_(n_state_ * n_particles_total_),
    state_next_(n_state_ * n_particles_total_),
    shared_(shared),
    internal_(internal),
    time_(time),
    dt_(dt),
    zero_every_(zero_every_vec<T>(shared_)),
    rng_(n_particles_total_, seed, deterministic) {
    // We don't check that the size is the same across all states;
    // this should be done by the caller (similarly, we don't check
    // that shared and internal have the same size).
  }

  auto run_to_time(real_type time) {
    const size_t n_steps = std::round(std::max(0.0, time - time_) / dt_);
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
                     zero_every_[i],
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
        *output = T::compare_data(time_, dt_, state_data + offset, *data,
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
  std::vector<zero_every_type<real_type>> zero_every_;
  mcstate::random::prng<rng_state_type> rng_;

  static void run_particle(real_type time, real_type dt, size_t n_steps,
                           const shared_state& shared, internal_state& internal,
                           const zero_every_type<real_type>& zero_every,
                           real_type * state, rng_state_type& rng_state,
                           real_type * state_next) {
    for (size_t i = 0; i < n_steps; ++i) {
      for (const auto& el : zero_every) {
        if (std::fmod(time, el.first) == 0) {
          for (const auto j : el.second) {
            state[j] = 0;
          }
        }
      }
      T::update(time + i * dt, dt, state, shared, internal, rng_state,
                state_next);
      std::swap(state, state_next);
    }
  }
};

}
