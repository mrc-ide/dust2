#pragma once

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
  adjoint_data(size_t n_state, size_t n_particles) :
    n_state_(n_state),
    n_adjoint_(0),
    n_particles_(n_particles),
    n_steps_(0) {
  }

  void init_history(size_t n_steps) {
    if (n_steps_ != n_steps) {
      state_.resize(n_state_ * n_particles_ * (n_steps + 1));
      n_steps_ = n_steps;
    }
  }

  void init_adjoint(size_t n_adjoint) {
    if (n_adjoint_ != n_adjoint) {
      adjoint_curr_.resize(n_adjoint * n_particles_);
      adjoint_next_.resize(adjoint_curr_.size());
      n_adjoint_ = n_adjoint;
    }
  }

  // TODO: it would be nicer to use iterators throughout I think,
  // though the end result will be identical?  Probably best to do
  // this throughout really and try and remove as many raw pointers as
  // possible.
  auto state() {
    return state_.data();
  }

  auto curr() {
    return adjoint_curr_.data();
  }

  auto next() {
    return adjoint_next_.data();
  }

  template <typename Iter>
  void gradient(Iter iter) {
    std::copy(adjoint_curr_.begin() + n_state_, adjoint_curr_.end(), iter);
  }

private:
  size_t n_state_;
  size_t n_adjoint_;
  size_t n_particles_;
  size_t n_steps_;
  std::vector<real_type> state_;
  std::vector<real_type> adjoint_curr_;
  std::vector<real_type> adjoint_next_;
};

}
