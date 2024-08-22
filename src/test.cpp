#include <dust2/filter_details.hpp>
#include <dust2/history.hpp>
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

// Simple driver for exercising the history saving outside of any
// particle filter.
[[cpp11::register]]
cpp11::sexp test_history_(cpp11::doubles r_time, cpp11::list r_state,
                          cpp11::sexp r_order,
                          cpp11::sexp r_index_group,
			  cpp11::sexp r_index_particle,
                          bool reorder) {
  const size_t n_times = r_time.size();
  cpp11::sexp el0 = r_state[0];
  auto r_dim = cpp11::as_cpp<cpp11::integers>(el0.attr("dim"));

  std::vector<size_t> all_groups;
  const size_t n_groups_all = r_dim[2];
  for (size_t i = 0; i < n_groups_all; ++i) {
    all_groups.push_back(i);
  }

  const size_t n_state = r_dim[0];
  const size_t n_particles = r_dim[1];

  const auto index_group = r_index_group == R_NilValue ? all_groups :
    dust2::r::check_index(r_index_group, all_groups.size(), "index_group");
  const size_t n_groups = index_group.size();

  std::vector<size_t> index_particle;
  const bool use_index_particle = r_index_particle != R_NilValue;
  if (use_index_particle) {
    dust2::r::check_length(r_index_particle, n_groups, "index_particle");
    index_particle = dust2::r::check_index(r_index_particle, n_particles,
					   "index_particle");
  }


  const size_t n_particles_out = index_particle.empty() ? n_particles : 1;

  dust2::history<double> h(n_state, n_particles, n_groups_all, n_times);
  for (size_t i = 0; i < static_cast<size_t>(r_state.size()); ++i) {
    if (r_order == R_NilValue) {
      h.add(r_time[i], REAL(r_state[i]), index_group);
    } else {
      cpp11::sexp el = cpp11::as_cpp<cpp11::list>(r_order)[i];
      if (el == R_NilValue) {
        h.add(r_time[i], REAL(r_state[i]), index_group);
      } else {
        h.add(r_time[i], REAL(r_state[i]), INTEGER(el), index_group);
      }
    }
  }

  const auto n_times_out = h.n_times();
  cpp11::writable::doubles ret_time(static_cast<int>(n_times_out));
  const size_t len = n_state * n_particles_out * n_groups * n_times_out;
  cpp11::writable::doubles ret_state(static_cast<int>(len));
  h.export_time(REAL(ret_time));
  h.export_state(REAL(ret_state), reorder, index_group, index_particle);

  if (use_index_particle) {
    dust2::r::set_array_dims(ret_state, {n_state, n_particles_out * n_groups, n_times_out});
  } else {
    dust2::r::set_array_dims(ret_state, {n_state, n_particles_out, n_groups, n_times_out});
  }

  return cpp11::writable::list{ret_time, ret_state};
}
