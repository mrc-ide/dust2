#pragma once

#include <vector>

template <typename real_type>
class adjoint_data {

public:
  adjoint_data(size_t n_state, size_t n_particles, size_t n_steps) :
    n_state_(n_state),
    n_particles_(n_particles),
    n_steps_(n_steps) {
  }

  void init_history() {
    if (state_.size() == 0) {
      state_.resize(n_state_ * n_particles_ * (n_steps_ + 1));
    }
  }

  void init_adjoint(size_t n_adjoint) {
    if (n_adjoint_ != n_adjoint) {
      adjoint_curr_.resize((n_state_ + n_adjoint) * n_particles_);
      adjoint_next_.resize(adjoint_curr_.size());
      n_adjoint_ = n_adjoint;
    }
  }

  // TODO: it would be nicer to use iterators throughout I think,
  // though the end result will be identical?  Probably best to do
  // this throughout really.
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
