#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

namespace dust2 {

// This class looks after all the bookkeeping for running an adjoint.
// It is unconditionally included within an unfilter, but because it
// can be quite large the memory is not allocated until it is actually
// used the first time.  There are two points of initialisation -
// setting up the history memory (which is large as it is all states
// by all times) and then the adjoint memory (which is quite small but
// requires the model to actually support adjoint running).
template <typename real_type>
class adjoint_data {
public:
  adjoint_data(size_t n_state, size_t n_particles, size_t n_groups) :
    n_state_(n_state),
    n_adjoint_(0),
    n_particles_(n_particles),
    n_groups_(n_groups),
    n_particles_total_(n_particles_ * n_groups_),
    n_steps_(0) {
  }

  void init_history(size_t n_steps) {
    if (n_steps_ != n_steps) {
      state_.resize(n_state_ * n_particles_total_ * (n_steps + 1));
      n_steps_ = n_steps;
    }
    reset();
  }

  void init_adjoint(size_t n_adjoint) {
    if (n_adjoint_ != n_adjoint) {
      adjoint_curr_.resize(n_adjoint * n_particles_total_);
      adjoint_next_.resize(adjoint_curr_.size());
      n_adjoint_ = n_adjoint;
    }
  }

  auto state() {
    return state_.data();
  }

  auto curr() {
    return adjoint_curr_.data();
  }

  auto next() {
    return adjoint_next_.data();
  }

  void reset() {
    std::fill(state_.begin(), state_.end(), 0);
    std::fill(adjoint_curr_.begin(), adjoint_curr_.end(), 0);
    std::fill(adjoint_next_.begin(), adjoint_next_.end(), 0);
  }

  template <typename Iter>
  void gradient(Iter iter, const std::vector<size_t>& index_group) {
    auto iter_src = adjoint_curr_.begin() + n_state_;
    const auto n_gradient = n_adjoint_ - n_state_;
    for (auto i : index_group) {
      for (size_t j = 0; j < n_particles_; ++j) {
        iter = std::copy_n(iter_src + n_adjoint_ * n_particles_ * i, n_gradient * n_particles_, iter);
      }
    }
  }

private:
  size_t n_state_;
  size_t n_adjoint_;
  size_t n_particles_;
  size_t n_groups_;
  size_t n_particles_total_;
  size_t n_steps_;
  std::vector<real_type> state_;
  std::vector<real_type> adjoint_curr_;
  std::vector<real_type> adjoint_next_;
};

}
