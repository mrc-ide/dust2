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
