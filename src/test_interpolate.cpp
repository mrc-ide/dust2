#include <cpp11.hpp>
#include <dust2/interpolate/interpolate.hpp>

[[cpp11::register]]
int test_interpolate_search(double target, std::vector<double> x) {
  return dust2::interpolate::internal::interpolate_search(target, x, false);
}

[[cpp11::register]]
double test_interpolate_constant1(std::vector<double> t, std::vector<double> y,
                                  double z) {
  auto obj = dust2::interpolate::InterpolateConstant<double>(t, y, "t", "y");
  return obj.eval(z);
}

[[cpp11::register]]
double test_interpolate_linear1(std::vector<double> t, std::vector<double> y,
                                  double z) {
  auto obj = dust2::interpolate::InterpolateLinear<double>(t, y, "t", "y");
  return obj.eval(z);
}

[[cpp11::register]]
double test_interpolate_spline1(std::vector<double> t, std::vector<double> y,
                                  double z) {
  auto obj = dust2::interpolate::InterpolateSpline<double>(t, y, "t", "y");
  return obj.eval(z);
}

[[cpp11::register]]
std::vector<double> test_interpolate_constant2(std::vector<double> t,
                                               cpp11::doubles r_y,
                                               double z) {
  std::vector<double> y(r_y.begin(), r_y.end());
  auto r_dim = cpp11::as_cpp<cpp11::integers>(r_y.attr("dim"));
  const dust2::array::dimensions<2> dim{static_cast<size_t>(r_dim[0]),
                                        static_cast<size_t>(r_dim[1])};
  std::vector<double> ret(dim.size);
  auto obj =
    dust2::interpolate::InterpolateConstantArray<double, 2>(t, y, dim, "t", "y");
  obj.eval(z, ret);
  return ret;
}

[[cpp11::register]]
std::vector<double> test_interpolate_linear2(std::vector<double> t,
                                               cpp11::doubles r_y,
                                               double z) {
  std::vector<double> y(r_y.begin(), r_y.end());
  auto r_dim = cpp11::as_cpp<cpp11::integers>(r_y.attr("dim"));
  const dust2::array::dimensions<2> dim{static_cast<size_t>(r_dim[0]),
                                        static_cast<size_t>(r_dim[1])};
  std::vector<double> ret(dim.size);
  auto obj =
    dust2::interpolate::InterpolateLinearArray<double, 2>(t, y, dim, "t", "y");
  obj.eval(z, ret);
  return ret;
}

[[cpp11::register]]
std::vector<double> test_interpolate_spline2(std::vector<double> t,
                                               cpp11::doubles r_y,
                                               double z) {
  std::vector<double> y(r_y.begin(), r_y.end());
  auto r_dim = cpp11::as_cpp<cpp11::integers>(r_y.attr("dim"));
  const dust2::array::dimensions<2> dim{static_cast<size_t>(r_dim[0]),
                                        static_cast<size_t>(r_dim[1])};
  std::vector<double> ret(dim.size);
  auto obj =
    dust2::interpolate::InterpolateSplineArray<double, 2>(t, y, dim, "t", "y");
  obj.eval(z, ret);
  return ret;
}
