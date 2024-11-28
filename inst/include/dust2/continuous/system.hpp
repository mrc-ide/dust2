#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <vector>
#include <dust2/continuous/control.hpp>
#include <dust2/continuous/delays.hpp>
#include <dust2/continuous/events.hpp>
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
    n_state_output_(do_n_state_output<T>(packing_state_)),
    n_state_ode_(n_state_ - n_state_output_),
    n_particles_(n_particles),
    n_groups_(shared.size()),
    n_threads_(n_threads),
    n_particles_total_(n_particles_ * n_groups_),

    dt_(dt),
    control_(control),

    state_(n_state_ * n_particles_total_),
    ode_internals_(n_particles_total_,
                   {n_state_ode_, control_.save_history || has_delays_}),

    // For reordering to work:
    state_other_(n_state_ * n_particles_total_),
    ode_internals_other_(ode_internals_),

    shared_(shared),
    internal_(internal),
    all_groups_(tools::integer_sequence(n_groups_)),

    time_(time),
    zero_every_(zero_every_vec<T>(shared_)),
    errors_(n_particles_total_),
    rng_(n_particles_total_, seed, deterministic),
    delays_(do_delays<T>(shared_)),
    events_(do_events<T>(shared_, internal_)),
    solver_(n_groups_ * n_threads_, {n_state_ode_, control_}),
    output_is_current_(n_groups_),
    requires_initialise_(n_groups_, true) {
    initialise_delays_();
  }

  template <typename mixed_time = typename dust2::properties<T>::is_mixed_time>
  typename std::enable_if<!mixed_time::value, void>::type
  run_to_time(real_type time, const std::vector<size_t>& index_group) {
    initialise_solver_(index_group);
    real_type * state_data = state_.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (auto group : index_group) {
      for (size_t particle = 0; particle < n_particles_; ++particle) {
        const auto thread = tools::thread_index();
        const auto i = thread * n_groups_ + group;
        const auto k = n_particles_ * group + particle;
        const auto offset = k * n_state_;
        real_type * y = state_data + offset;
        try {
          solver_[i].run(time_, time, y, zero_every_[group], events_[i],
                         ode_internals_[k],
                         rhs_(particle, group, thread));
        } catch (std::exception const& e) {
          errors_.capture(e, k);
        }
      }
    }
    errors_.report();
    time_ = time;
    update_output_is_current(index_group, false);
  }

  template <typename mixed_time = typename dust2::properties<T>::is_mixed_time>
  typename std::enable_if<mixed_time::value, void>::type
  run_to_time(real_type time, const std::vector<size_t>& index_group) {
    initialise_solver_(index_group);
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
    for (auto group : index_group) {
      for (size_t particle = 0; particle < n_particles_; ++particle) {
        const auto thread = tools::thread_index();
        const auto i = thread * n_groups_ + group;
        const auto k = n_particles_ * group + particle;
        const auto offset = k * n_state_;
        auto& rng_state = rng_.state(k);
        real_type * y = state_data + offset;
        real_type * y_other = state_other_data + offset;
        try {
          real_type t1 = time_;
          for (size_t step = 0; step < n_steps; ++step) {
            const real_type t0 = t1;
            t1 = (step == n_steps - 1) ? time : time_ + step * dt_;
            solver_[i].run(t0, t1, y, zero_every_[group], events_[i],
                           ode_internals_[k],
                           rhs_(particle, group, thread));
            std::copy_n(y, n_state_ode_, y_other);
            update_(particle, group, thread,
                    time, y, rng_state, y_other);
            std::swap(y, y_other);
            solver_[i].initialise(time_, y, ode_internals_[k],
                                  rhs_(particle, group, thread));
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
    update_output_is_current(index_group, false);
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
        } catch (std::exception const& e) {
          errors_.capture(e, k);
        }
      }
    }
    errors_.report();
    mark_requires_initialise(index_group);
    // Assume not current, because most models would want to call output here()
    update_output_is_current(index_group, false);
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
    mark_requires_initialise(index_group);
    // I'm not sure what is best here; this (and to a degree
    // T::initial) are the two places where we might end up with
    // inconsistent output (e.g., the user has set a state that
    // includes output but the output is wrong).  This approach gives
    // them some flexibility at least.
    bool update_is_current = tools::is_trivial_index(index_state, n_state_);
    update_output_is_current(index_group, update_is_current);
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
        const auto n_state_copy =
          output_is_current_[i] ? n_state_ : n_state_ode_;
        std::copy_n(state_.begin() + k_from * n_state_,
                    n_state_copy,
                    state_other_.begin() + k_to * n_state_);
        ode_internals_other_[k_to] = ode_internals_[k_from];
      }
    }
    std::swap(state_, state_other_);
    std::swap(ode_internals_, ode_internals_other_);
  }

  const auto& state() {
    update_output();
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
    mark_requires_initialise(all_groups_);
    // We should set output_is_current here but I will wait until
    // updating time_ to make it vary by group. Practically the next
    // thing anyone does after this is to update initial conditions so
    // this is fine for now.
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
    // TODO: we might want to update the delays here, too.  Always
    // doing so would be the safest, but we might just forbid it and
    // document that as such (mrc-5993)
    requires_initialise_[i] = true;
  }

  template <typename IterData, typename IterOutput>
  void compare_data(IterData data,
                    const size_t n_groups_data,
                    const std::vector<size_t>& index_group,
                    IterOutput output) {
    compare_data(data, n_groups_data, state_.data(), index_group, output);
  }

  template <typename IterData, typename IterOutput>
  void compare_data(IterData data,
                    const size_t n_groups_data,
                    const real_type * state,
                    const std::vector<size_t>& index_group,
                    IterOutput output) {
    bool data_is_shared = n_groups_data == 1; // see discrete/system.hpp
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (auto i : index_group) {
      for (size_t j = 0; j < n_particles_; ++j) {
        const auto data_i = data + (data_is_shared ? 0 : i);
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
  size_t n_state_output_;
  size_t n_state_ode_;
  size_t n_particles_;
  size_t n_groups_;
  size_t n_threads_;
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
  std::vector<ode::delay_result_type<real_type>> delay_result_;
  std::vector<size_t> all_groups_;
  real_type time_;
  std::vector<zero_every_type<real_type>> zero_every_;
  dust2::utils::errors errors_;
  monty::random::prng<rng_state_type> rng_;
  std::vector<ode::delays<real_type>> delays_;
  std::vector<ode::events_type<real_type>> events_;
  std::vector<ode::solver<real_type>> solver_;
  std::vector<bool> output_is_current_;
  std::vector<bool> requires_initialise_;
  static constexpr bool has_delays_ = properties<T>::has_delays::value;
  static constexpr bool rhs_uses_delays_ = properties<T>::rhs_uses_delays;
  static constexpr bool output_uses_delays_ = properties<T>::output_uses_delays;

  auto rhs_(size_t particle, size_t group, size_t thread) {
    const size_t i = thread * n_groups_ + group;
    const size_t j = group * n_particles_ + particle;
    const auto& shared = shared_[group];
    auto& internal = internal_[i];
    const auto& history = ode_internals_[j].history_values;
    auto& delays = delays_[group];
    auto& delay_result = delay_result_[i];
    return [&](real_type t, const real_type* y, real_type* dydt) {
      if constexpr (rhs_uses_delays_) {
        delays.eval(t, history, delay_result);
        T::rhs(t, y, shared, internal, delay_result, dydt);
      } else {
        T::rhs(t, y, shared, internal, dydt);
      }
    };
  }

  void output_(size_t particle, size_t group, size_t thread) {
    if (!output_is_current_[group]) {
      const size_t i = thread * n_groups_ + group;
      const size_t j = group * n_particles_ + particle;
      const auto& shared = shared_[group];
      auto& internal = internal_[i];
      auto& delays = delays_[group];
      const auto& history = ode_internals_[j].history_values;

      real_type * y = state_.data() + j * n_state_;
      if constexpr (output_uses_delays_) {
        delays.eval(time_, history, delay_result_[i]);
        T::output(time_, y, shared, internal, delay_result_[i]);
      } else {
        T::output(time_, y, shared, internal);
      }
    }
  }

  void update_(size_t particle, size_t group, size_t thread,
               real_type time, const double * y, rng_state_type rng_state,
               double * y_other) {
    const size_t i = thread * n_groups_ + group;
    const auto& shared = shared_[group];
    auto& internal = internal_[i];

    T::update(time, dt_, y, shared, internal, rng_state, y_other);
  }

  void initialise_delays_() {
    if constexpr (has_delays_) {
      // Create the space we need to save results - like internal we
      // do this once per thread as we don't retain contents across
      // time steps.
      delay_result_.reserve(n_particles_ * n_groups_);
      for (size_t i = 0, k = 0; i < n_threads_; ++i) {
        for (size_t group = 0; group < n_groups_; ++group, ++k) {
          delay_result_.push_back(delays_[group].result());
          solver_[k].control().step_size_max =
            delays_[group].step_size_max(control_.step_size_max,
                                         rhs_uses_delays_);
        }
      }
      for (size_t i = 0, k = 0; i < n_groups_; ++i) {
        for (size_t j = 0; j < n_particles_; ++j, ++k) {
          ode_internals_[k].history_values.set_index(delays_[i].index());
          ode_internals_other_[k].history_values.set_index(delays_[i].index());
        }
      }
      control_.save_history = true;
    } else {
      delay_result_.resize(n_particles_ * n_groups_);
    }
  }

  void initialise_solver_(std::vector<size_t> index_group) {
    bool requires_initialisation = false;
    for (auto i : index_group) {
      if (requires_initialise_[i]) {
        requires_initialisation = true;
        break;
      }
    }
    if (!requires_initialisation) {
      return;
    }
    errors_.reset();
    real_type * state_data = state_.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (auto group : index_group) {
      for (size_t particle = 0; particle < n_particles_; ++particle) {
        if (requires_initialise_[group]) {
          const auto thread = tools::thread_index();
          const auto i = thread * n_groups_ + group;
          const auto k = n_particles_ * group + particle;
          const auto offset = k * n_state_;
          real_type * y = state_data + offset;
          try {
            solver_[i].initialise(time_, y, ode_internals_[k],
                                  rhs_(particle, group, thread));
          } catch (std::exception const& e) {
            errors_.capture(e, k);
          }
        }
      }
    }
    errors_.report();
    for (auto i : index_group) {
      requires_initialise_[i] = false;
    }
  }

  // Default implementation does nothing
  template <typename has_output = typename dust2::properties<T>::has_output>
  typename std::enable_if<!has_output::value, void>::type
  update_output() {
  }

  template <typename has_output = typename dust2::properties<T>::has_output>
  typename std::enable_if<has_output::value, void>::type
  update_output() {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) collapse(2)
#endif
    for (size_t group = 0; group < n_groups_; ++group) {
      for (size_t particle = 0; particle < n_particles_; ++particle) {
        const auto thread = tools::thread_index();
        const auto k = n_particles_ * group + particle;
        try {
          output_(particle, group, thread);
        } catch (std::exception const& e) {
          errors_.capture(e, k);
        }
      }
    }
    errors_.report();
    update_output_is_current({}, true);
  }

  void update_output_is_current(std::vector<size_t> index_group, bool value) {
    if (index_group.empty()) {
      std::fill(output_is_current_.begin(), output_is_current_.end(), value);
    } else {
      for (auto i : index_group) {
        output_is_current_[i] = value;
      }
    }
  }

  void mark_requires_initialise(std::vector<size_t> index_group) {
    if (index_group.empty()) {
      std::fill(requires_initialise_.begin(), requires_initialise_.end(), true);
    } else {
      for (auto i : index_group) {
        requires_initialise_[i] = true;
      }
    }
  }
};

}
