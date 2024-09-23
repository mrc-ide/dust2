#include <dust2/array.hpp>
#include <cpp11.hpp>

[[cpp11::register]]
double test_sum(cpp11::sexp r_x, cpp11::sexp r_m) {
  const double* x = REAL(r_x);
  cpp11::integers r_dim = r_x.attr("dim") == R_NilValue ?
    cpp11::writable::integers{LENGTH(r_x)} :
    cpp11::as_cpp<cpp11::integers>(r_x.attr("dim"));
  const size_t rank = LENGTH(r_dim);
  double ret = 0;

  std::vector<size_t> m;
  if (r_m != R_NilValue) {
    auto tmp = cpp11::as_cpp<cpp11::integers>(r_m);
    for (int i = 0; i < tmp.size(); ++i) {
      m.push_back(tmp[i] - 1);
    }
  }

  if (rank == 1) {
    dust2::array::dimensions<1> dim{static_cast<size_t>(r_dim[0])};
    if (r_m == R_NilValue) {
      ret = dust2::array::sum<double>(x, dim);
    } else {
      ret = dust2::array::sum<double>(x, dim, {m[0], m[1]});
    }
  } else if (rank == 2) {
    dust2::array::dimensions<2> dim{static_cast<size_t>(r_dim[0]),
                                    static_cast<size_t>(r_dim[1])};
    if (r_m == R_NilValue) {
      ret = dust2::array::sum<double>(x, dim);
    } else {
      ret = dust2::array::sum<double>(x, dim, {m[0], m[1]}, {m[2], m[3]});
    }
  } else if (rank == 3) {
    dust2::array::dimensions<3> dim{static_cast<size_t>(r_dim[0]),
                                    static_cast<size_t>(r_dim[1]),
                                    static_cast<size_t>(r_dim[2])};
    if (r_m == R_NilValue) {
      ret = dust2::array::sum<double>(x, dim);
    } else {
      ret = dust2::array::sum<double>(x, dim, {m[0], m[1]}, {m[2], m[3]}, {m[4], m[5]});
    }
  }

  return ret;
}


[[cpp11::register]]
double test_min(cpp11::sexp r_x, cpp11::sexp r_m) {
  const double* x = REAL(r_x);
  cpp11::integers r_dim = r_x.attr("dim") == R_NilValue ?
    cpp11::writable::integers{LENGTH(r_x)} :
    cpp11::as_cpp<cpp11::integers>(r_x.attr("dim"));
  const size_t rank = LENGTH(r_dim);
  double ret = 0;

  std::vector<size_t> m;
  if (r_m != R_NilValue) {
    auto tmp = cpp11::as_cpp<cpp11::integers>(r_m);
    for (int i = 0; i < tmp.size(); ++i) {
      m.push_back(tmp[i] - 1);
    }
  }

  if (rank == 1) {
    dust2::array::dimensions<1> dim{static_cast<size_t>(r_dim[0])};
    if (r_m == R_NilValue) {
      ret = dust2::array::min<double>(x, dim);
    } else {
      ret = dust2::array::min<double>(x, dim, {m[0], m[1]});
    }
  } else if (rank == 2) {
    dust2::array::dimensions<2> dim{static_cast<size_t>(r_dim[0]),
                                    static_cast<size_t>(r_dim[1])};
    if (r_m == R_NilValue) {
      ret = dust2::array::min<double>(x, dim);
    } else {
      ret = dust2::array::min<double>(x, dim, {m[0], m[1]}, {m[2], m[3]});
    }
  } else if (rank == 3) {
    dust2::array::dimensions<3> dim{static_cast<size_t>(r_dim[0]),
                                    static_cast<size_t>(r_dim[1]),
                                    static_cast<size_t>(r_dim[2])};
    if (r_m == R_NilValue) {
      ret = dust2::array::min<double>(x, dim);
    } else {
      ret = dust2::array::min<double>(x, dim, {m[0], m[1]}, {m[2], m[3]}, {m[4], m[5]});
    }
  }

  return ret;
}
