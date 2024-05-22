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

// Simple driver for exercising the history saving outside of any
// particle filter.
[[cpp11::register]]
cpp11::sexp test_history(cpp11::doubles r_time, cpp11::list r_state,
                         cpp11::sexp r_order, bool reorder) {
  const size_t n_times = r_time.size();
  cpp11::sexp el0 = r_state[0];

  auto r_dim = cpp11::as_cpp<cpp11::integers>(el0.attr("dim"));
  const size_t n_state = r_dim[0];
  const size_t n_particles = r_dim[1];
  const size_t n_groups = r_dim[2];

  dust2::history<double> h(n_state, n_particles, n_groups, n_times);
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

  cpp11::writable::doubles ret_time(static_cast<int>(h.size_time()));
  cpp11::writable::doubles ret_state(static_cast<int>(h.size_state()));
  h.export_time(REAL(ret_time));
  h.export_state(REAL(ret_state), reorder);
  dust2::r::set_array_dims(ret_state, {n_state, n_particles, n_groups, h.size_time()});

  return cpp11::writable::list{ret_time, ret_state};
}
