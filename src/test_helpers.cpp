#include <dust2/r/helpers.hpp>

[[cpp11::register]]
bool test_check_dimensions(cpp11::sexp value, cpp11::integers r_dim,
                           const char * name) {
  const size_t rank = r_dim.size();
  if (rank == 1) {
    const dust2::array::dimensions<1> dim{static_cast<size_t>(r_dim[0])};
    dust2::r::check_dimensions(value, dim, name);
  } else if (rank == 2) {
    const dust2::array::dimensions<2> dim{static_cast<size_t>(r_dim[0]),
                                          static_cast<size_t>(r_dim[1])};
    dust2::r::check_dimensions(value, dim, name);
  } else if (rank == 3) {
    const dust2::array::dimensions<3> dim{static_cast<size_t>(r_dim[0]),
                                          static_cast<size_t>(r_dim[1]),
                                          static_cast<size_t>(r_dim[2])};
    dust2::r::check_dimensions(value, dim, name);
  }
  return true;
}


[[cpp11::register]]
std::vector<size_t> test_read_dimensions(cpp11::list value, int rank, const char * name) {
  if (rank == 1) {
    const auto arr = dust2::r::read_dimensions<1>(value, name).dim;
    return std::vector<size_t>(arr.begin(), arr.end());
  } else if (rank == 2) {
    const auto arr = dust2::r::read_dimensions<2>(value, name).dim;
    return std::vector<size_t>(arr.begin(), arr.end());
  } else {
    const auto arr = dust2::r::read_dimensions<3>(value, name).dim;
    return std::vector<size_t>(arr.begin(), arr.end());
  }
}
