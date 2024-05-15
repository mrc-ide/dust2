#include <dust2/filter_details.hpp>

#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>

[[cpp11::register]]
cpp11::integers test_resample_weight(std::vector<double> w, double u) {
  const size_t n_particles = w.size();
  std::vector<size_t> idx(n_particles);
  dust2::details::resample_weight(n_particles, w.begin(), u, idx.begin());
  cpp11::writable::integers ret(idx.begin(), idx.end());
  return ret;
}
