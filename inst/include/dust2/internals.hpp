#pragma once

#include <algorithm>
#include <vector>

namespace dust2 {
namespace internals {

template <typename real_type, typename Iter>
void set_state(std::vector<real_type>& state,
               Iter iter,
               const size_t n_state,
               const size_t n_particles,
               const size_t n_groups,
               const std::vector<size_t>& index_state,
               const std::vector<size_t>& index_particle,
               const std::vector<size_t>& index_group,
               const bool recycle_particle,
               const bool recycle_group,
               const size_t n_threads) {
  const bool use_index_state = !index_state.empty();
  const bool use_index_particle = !index_particle.empty();
  const bool use_index_group = !index_group.empty();

  bool do_simple_copy =
    !use_index_state && !use_index_particle && !use_index_group &&
    !recycle_particle && !recycle_group;
  if (do_simple_copy) {
    std::copy_n(iter, state.size(), state.begin());
  } else {
    const auto n_state_in =
      use_index_state ? index_state.size() : n_state;
    const auto n_particles_in =
      use_index_particle ? index_particle.size() : n_particles;
    const auto n_groups_in =
      use_index_group ? index_group.size() : n_groups;
    const auto offset_src_particle = recycle_particle ? 0 : n_state_in;
    const auto offset_src_group = recycle_group ? 0 :
      (n_state_in * (recycle_particle ? 1 : n_particles_in));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads) collapse(2)
#endif
    for (size_t i = 0; i < n_groups_in; ++i) {
      for (size_t j = 0; j < n_particles_in; ++j) {
        const auto offset_src =
          i * offset_src_group + j * offset_src_particle;
        const auto i_dst = use_index_group ? index_group[i] : i;
        const auto j_dst = use_index_particle ? index_particle[j] : j;
        const auto offset_dst = (n_particles * i_dst + j_dst) * n_state;
        const auto iter_src = iter + offset_src;
        auto iter_dst = state.begin() + offset_dst;
        if (use_index_state) {
          for (size_t k = 0; k < n_state_in; ++k) {
            *(iter_dst + index_state[k]) = *(iter_src + k);
          }
        } else {
          std::copy_n(iter_src, n_state, iter_dst);
        }
      }
    }
  }
}

}
}
