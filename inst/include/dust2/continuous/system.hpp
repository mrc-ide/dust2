#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <vector>
#include <dust2/continuous/control.hpp>
#include <dust2/continuous/solver.hpp>
#include <dust2/errors.hpp>
#include <dust2/internals.hpp>
#include <dust2/packing.hpp>
#include <dust2/properties.hpp>
#include <dust2/tools.hpp>
#include <dust2/zero.hpp>
#include <monty/random/random.hpp>

namespace dust2 {

template <typename T>
class dust_continuous {
public:
  using system_type = T;
  using real_type = typename T::real_type;
  using rng_state_type = typename T::rng_state_type;
  using shared_state = typename T::shared_state;
  using internal_state = typename T::internal_state;
  using rng_int_type = typename rng_state_type::int_type;

  dust_continuous(std::vector<shared_state> shared,
                  std::vector<internal_state> internal,
                  real_type time,
                  real_type dt,
                  const ode::control<real_type> control,
                  size_t n_particles, // per group
                  const std::vector<rng_int_type>& seed,
                  bool deterministic, size_t n_threads) :
    packing_state_(T::packing_state(shared[0])),
    packing_gradient_(do_packing_gradient<T>(shared[0])),
    n_state_(packing_state_.size()),
    n_particles_(n_particles),
    n_groups_(shared.size()),
    n_particles_total_(n_particles_ * n_groups_),

    dt_(dt),
    control_(control),

    state_(n_state_ * n_particles_total_),
    ode_internals_(n_particles_total_, n_state_),

    // For reordering to work:
    state_other_(n_state_ * n_particles_total_),
    ode_internals_other_(n_particles_total_, n_state_),

    shared_(shared),
    internal_(internal),
    all_groups_(n_groups_),

    time_(time),
    zero_every_(zero_every_vec<T>(shared_)),
    errors_(n_particles_total_),
    rng_(n_particles_total_, seed, deterministic),
    solver_(n_state_, control_),
    n_threads_(n_threads) {
    // TODO: above, filter rng states need adding here too, or
    // somewhere at least (we might move the filter elsewhere though,
    // in which case that particular bit of weirdness goes away).

    // We don't check that the size is the same across all states;
    // this should be done by the caller (similarly, we don't check
    // that shared and internal have the same size).
    for (size_t i = 0; i < n_groups_; ++i) {
      all_groups_[i] = i;
    }
  }

  template <typename mixed_time = typename dust2::properties<T>::is_mixed_time>
  typename std::enable_if<!mixed_time::value, void>::type
  run_to_time(real_type time, const std::vector<size_t>& index_group) {
    real_type * state_data = state_.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (auto i : index_group) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        auto& internal_i = internal_[tools::thread_index() * n_groups_ + i];
        real_type * y = state_data + offset;
        try {
          solver_.run(time_, time, y, zero_every_[i],
                      ode_internals_[k],
                      rhs_(shared_[i], internal_i));
        } catch (std::exception const& e) {
          errors_.capture(e, k);
        }
      }
    }
    errors_.report();
    time_ = time;
  }

  template <typename mixed_time = typename dust2::properties<T>::is_mixed_time>
  typename std::enable_if<mixed_time::value, void>::type
  run_to_time(real_type time, const std::vector<size_t>& index_group) {
    if (dt_ == 0) {
      run_to_time<std::false_type>(time, index_group);
      return;
    }
    real_type * state_data = state_.data();
    real_type * state_other_data = state_other_.data();
    const size_t n_steps = (time - time_) / dt_;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (auto i : index_group) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        auto& internal_i = internal_[tools::thread_index() * n_groups_ + i];
        auto& rng_state = rng_.state(k);
        real_type * y = state_data + offset;
        real_type * y_other = state_other_data + offset;
        try {
          real_type t1 = time_;
          for (size_t step = 0; step < n_steps; ++step) {
            const real_type t0 = t1;
            t1 = (step == n_steps - 1) ? time : time_ + step * dt_;
            solver_.run(t0, t1, y, zero_every_[i],
                        ode_internals_[k],
                        rhs_(shared_[i], internal_i));
            std::copy_n(y, n_state_, y_other);
            T::update(t1, dt_, y, shared_[i], internal_i, rng_state, y_other);
            std::swap(y, y_other);
            solver_.initialise(time_, y, ode_internals_[k],
                               rhs_(shared_[i], internal_i));
          }
          } catch (std::exception const& e) {
          errors_.capture(e, k);
        }
      }
    }
    errors_.report();
    if (n_steps % 2 == 1) {
      std::swap(state_, state_other_);
    }
    time_ = time;
  }

  void run_to_time(real_type time,
                   const std::vector<size_t>& index_group,
                   real_type *state_history) {
    throw std::runtime_error("Write run_to_time() with saved history");
  }

  void set_state_initial(const std::vector<size_t>& index_group) {
    errors_.reset();
    real_type * state_data = state_.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (auto i : index_group) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        auto& internal_i = internal_[tools::thread_index() * n_groups_ + i];
        real_type * y = state_data + offset;
        try {
          T::initial(time_, shared_[i], internal_i,
                     rng_.state(k), y);
          solver_.initialise(time_, y, ode_internals_[k],
                             rhs_(shared_[i], internal_i));
        } catch (std::exception const& e) {
          errors_.capture(e, k);
        }
      }
    }
    errors_.report();
  }

  template <typename Iter>
  void set_state(Iter iter,
                 const std::vector<size_t>& index_state,
                 const std::vector<size_t>& index_particle,
                 const std::vector<size_t>& index_group,
                 bool recycle_particle,
                 bool recycle_group) {
    errors_.reset();
    dust2::internals::set_state(state_, iter,
                                n_state_, n_particles_, n_groups_,
                                index_state, index_particle, index_group,
                                recycle_particle, recycle_group,
                                n_threads_);
    initialise_solver_(index_group.empty() ? all_groups_ : index_group);
  }

  // iter here is an iterator to our *reordering index*, which will be
  // in terms of particles.  This probably warrants a better name.
  template <typename Iter>
  void reorder(Iter iter, const std::vector<size_t>& index_group) {
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
    state_other_.resize(state_.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (auto i : index_group) {
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

  auto& all_groups() const {
    return all_groups_;
  }

  auto& packing_state() const {
    return packing_state_;
  }

  auto& packing_gradient() const {
    return packing_gradient_;
  }

  auto& shared() const {
    return shared_;
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
  void compare_data(IterData data,
                    const std::vector<size_t>& index_group,
                    IterOutput output) {
    compare_data(data, state_.data(), index_group, output);
  }

  template <typename IterData, typename IterOutput>
  void compare_data(IterData data, const real_type * state,
                    const std::vector<size_t>& index_group,
                    IterOutput output) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (auto i : index_group) {
      for (size_t j = 0; j < n_particles_; ++j) {
        auto data_i = data + i;
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        auto& internal_i = internal_[tools::thread_index() * n_groups_ + i];
        auto output_ij = output + k;
        try {
          *output_ij = T::compare_data(time_, state + offset, *data_i,
                                       shared_[i], internal_i, rng_.state(k));
        } catch (std::exception const& e) {
          errors_.capture(e, k);
        }
      }
    }
    errors_.report();
  }

  // This is just used for debugging
  const auto& ode_internals() const {
    return ode_internals_;
  }

  bool errors_pending() const {
    return errors_.unresolved();
  }

private:
  dust2::packing packing_state_;
  dust2::packing packing_gradient_;
  size_t n_state_;
  size_t n_particles_;
  size_t n_groups_;
  size_t n_particles_total_;
  real_type dt_;
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
  std::vector<size_t> all_groups_;
  real_type time_;
  std::vector<zero_every_type<real_type>> zero_every_;
  dust2::utils::errors errors_;
  monty::random::prng<rng_state_type> rng_;
  ode::solver<real_type> solver_;
  size_t n_threads_;

  static auto rhs_(const shared_state& shared, internal_state& internal) {
    return [&](real_type t, const real_type* y, real_type* dydt) {
             T::rhs(t, y, shared, internal, dydt);
           };
  }

  void initialise_solver_(std::vector<size_t> index_group) {
    errors_.reset();
    real_type * state_data = state_.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (auto i : index_group) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto k = n_particles_ * i + j;
        const auto offset = k * n_state_;
        auto& internal_i = internal_[tools::thread_index() * n_groups_ + i];
        real_type * y = state_data + offset;
        try {
          solver_.initialise(time_, y, ode_internals_[k],
                             rhs_(shared_[i], internal_i));
        } catch (std::exception const& e) {
          errors_.capture(e, k);
        }
      }
    }
    errors_.report();
  }
};

}
