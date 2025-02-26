#include <dust2/filter_details.hpp>
#include <dust2/trajectories.hpp>
#include <dust2/r/helpers.hpp>

#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/list.hpp>


[[cpp11::register]]
cpp11::integers test_resample_weight(std::vector<double> w, double u) {
  const size_t n_particles = w.size();
  std::vector<size_t> idx(n_particles);
  dust2::details::resample_weight(n_particles, w.begin(), u, idx.begin());
  cpp11::writable::integers ret(idx.begin(), idx.end());
  return ret;
}

[[cpp11::register]]
cpp11::list test_scale_log_weights(std::vector<double> w) {
  const auto res = dust2::details::scale_log_weights<double>(w.size(), w.begin());
  return cpp11::writable::list{cpp11::as_sexp(res), cpp11::as_sexp(w)};
}

// Simple driver for exercising the trajectory saving outside of any
// particle filter.
[[cpp11::register]]
cpp11::sexp test_trajectories_(cpp11::doubles r_time,
                               cpp11::list r_state,
                               cpp11::sexp r_order,
                               cpp11::sexp r_index_state,
                               cpp11::sexp r_index_group,
                               cpp11::sexp r_select_particle,
                               cpp11::sexp r_times_snapshot,
                               bool save_state,
                               bool reorder) {
  const size_t n_times = r_time.size();
  cpp11::sexp el0 = r_state[0];
  auto r_dim = cpp11::as_cpp<cpp11::integers>(el0.attr("dim"));
  const size_t n_state = r_dim[0];
  const size_t n_particles = r_dim[1];
  const size_t n_groups = r_dim[2];

  std::vector<size_t> all_groups;
  for (size_t i = 0; i < n_groups; ++i) {
    all_groups.push_back(i);
  }
  const auto index_state =
    dust2::r::check_index(r_index_state, n_state, "index_state");
  const auto index_group = r_index_group == R_NilValue ? all_groups :
    dust2::r::check_index(r_index_group, all_groups.size(), "index_group");

  std::vector<size_t> select_particle;
  const bool use_select_particle = r_select_particle != R_NilValue;
  if (use_select_particle) {
    dust2::r::check_length(r_select_particle, index_group.size(),
                           "select_particle");
    select_particle = dust2::r::check_index(r_select_particle, n_particles,
                                            "select_particle");
  }

  const auto times_snapshot = r_times_snapshot == R_NilValue ?
    std::vector<double>() :
    cpp11::as_cpp<std::vector<double>>(r_times_snapshot);

  dust2::trajectories<double> h(n_state, n_particles, n_groups, n_times);
  h.set_index_and_reset(index_state, index_group, save_state, times_snapshot);
  for (size_t i = 0; i < static_cast<size_t>(r_state.size()); ++i) {
    if (r_order == R_NilValue) {
      h.add(r_time[i], REAL(r_state[i]));
    } else {
      cpp11::sexp el = cpp11::as_cpp<cpp11::list>(r_order)[i];
      if (el == R_NilValue) {
        h.add(r_time[i], REAL(r_state[i]));
      } else {
        h.add(r_time[i], REAL(r_state[i]), INTEGER(el));
      }
    }
  }

  const size_t n_particles_out = select_particle.empty() ? n_particles : 1;
  const auto n_state_out = h.n_state();
  const auto n_groups_out = h.n_groups();
  const auto n_times_out = h.n_times();
  cpp11::writable::doubles ret_time(static_cast<int>(n_times_out));

  cpp11::sexp ret_state = R_NilValue;
  cpp11::sexp ret_snapshots = R_NilValue;

  if (save_state) {
    const size_t len = n_state_out * n_particles_out * n_groups_out * n_times_out;
    cpp11::writable::doubles arr(static_cast<int>(len));
    h.export_time(REAL(ret_time));
    h.export_state(REAL(arr), reorder, select_particle);

    if (use_select_particle) {
      dust2::r::set_array_dims(arr, {n_state_out, n_particles_out * n_groups_out, n_times_out});
    } else {
      dust2::r::set_array_dims(arr, {n_state_out, n_particles_out, n_groups_out, n_times_out});
    }
    ret_state = arr;
  }

  if (!times_snapshot.empty()) {
    const auto n_snapshots = h.n_snapshots();
    const size_t len = n_state * n_particles_out * n_groups_out * n_snapshots;
    cpp11::writable::doubles arr(static_cast<int>(len));
    h.export_snapshots(REAL(arr), reorder, select_particle);

    if (use_select_particle) {
      dust2::r::set_array_dims(arr, {n_state, n_particles_out * n_groups_out, n_snapshots});
    } else {
      dust2::r::set_array_dims(arr, {n_state, n_particles_out, n_groups_out, n_snapshots});
    }
    ret_snapshots = arr;
  }

  return cpp11::writable::list{ret_time, ret_state, ret_snapshots};
}
