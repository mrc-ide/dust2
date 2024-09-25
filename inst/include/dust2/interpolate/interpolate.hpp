#pragma once

#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <vector>

#include <dust2/array.hpp>
#include <dust2/interpolate/spline.hpp>

namespace dust2 {
namespace interpolate {

namespace internal {

// TODO: add a hint here, optionally
template <typename T>
size_t interpolate_search(T target, std::vector<T> container,
                          bool allow_extrapolate_rhs) {
  if (target < container[0]) {
    const auto delta = container[0] - target;
    std::stringstream msg;
    msg << "Tried to interpolate at time = " << target <<
      ", which is " << delta << " before the first time (" <<
      container[0] << ")";
    throw std::runtime_error(msg.str());
  }
  const auto n = container.size();
  auto lower = std::lower_bound(container.begin(), container.end(), target);

  if (lower == container.end()) {
    if (allow_extrapolate_rhs || container[n - 1] == target) {
      return n - 1;
    } else {
      std::stringstream msg;
      const auto delta = target - container[n - 1];
      msg << "Tried to interpolate at time = " << target <<
        ", which is " << delta << " after the last time (" <<
        container[n - 1] << ")";
      throw std::runtime_error(msg.str());
    }
  }

  const auto i = std::distance(container.begin(), lower);
  return container[i] != target ? i - 1 : i;
}

template <typename T>
void validate_time(const std::vector<T>& t, size_t len_y,
                   const char * name_t, const char * name_y) {
  if (t.size() != len_y) {
    std::stringstream msg;
    msg << "Time variable '" << name_t << "' and interpolation target '" <<
      name_y << "' must have the same length, but do not (" <<
      t.size() << " vs " << len_y << ")";
    throw std::runtime_error(msg.str());
  }
  for (size_t i = 1; i < t.size(); ++i) {
    if (t[i - 1] >= t[i]) {
      std::stringstream msg;
      msg << "Time variable '" << name_t << "' must be strictly increasing " <<
        "but was not at index " << i;
      throw std::runtime_error(msg.str());
    }
  }
}

}

template <typename T>
class InterpolateConstant {
private:
  std::vector<T> t_;
  std::vector<T> y_;
public:
  InterpolateConstant(const std::vector<T>& t, const std::vector<T>& y,
                      const char * name_t, const char * name_y) :
    t_(t), y_(y) {
    internal::validate_time(t_, y_.size(), name_t, name_y);
  }

  InterpolateConstant() {}

  T eval(T z) const {
    const auto i = internal::interpolate_search(z, t_, true);
    return y_[i];
  }
};

template <typename T, size_t rank>
class InterpolateConstantArray {
private:
  std::vector<T> t_;
  std::vector<T> y_;
  size_t n_;
public:
  InterpolateConstantArray(const std::vector<T> t, const std::vector<T>& y,
                           const array::dimensions<rank>& dim,
                           const char * name_t, const char * name_y) :
    t_(t), y_(y), n_(dim.size) {
    const size_t len_y = y_.size() / n_;
    internal::validate_time(t, len_y, name_t, name_y);
  }

  InterpolateConstantArray() {}

  template <typename Iterator>
  void eval(T z, Iterator dest) const {
    const auto i = internal::interpolate_search(z, t_, true);
    std::copy_n(y_.begin() + n_ * i, n_, dest);
  }
};

template <typename T>
class InterpolateLinear {
private:
  std::vector<T> t_;
  std::vector<T> y_;
public:
  InterpolateLinear(const std::vector<T>& t, const std::vector<T>& y,
                    const char * name_t, const char * name_y) :
    t_(t), y_(y) {
    internal::validate_time(t_, y_.size(), name_t, name_y);
  }

  InterpolateLinear() {}

  T eval(T z) const {
    auto i = internal::interpolate_search(z, t_, false);
    if (i == t_.size() - 1) {
      return y_[i];
    }
    T t0 = t_[i], t1 = t_[i + 1], y0 = y_[i], y1 = y_[i + 1];
    return y0 + (y1 - y0) * (z - t0) / (t1 - t0);
  }
};

template <typename T, size_t rank>
class InterpolateLinearArray {
private:
  std::vector<T> t_;
  std::vector<T> y_;
  size_t n_;
public:
  InterpolateLinearArray(const std::vector<T> t, const std::vector<T>& y,
                           const array::dimensions<rank>& dim,
                           const char * name_t, const char * name_y) :
    t_(t), y_(y), n_(dim.size) {
    const size_t len_y = y_.size() / n_;
    internal::validate_time(t, len_y, name_t, name_y);
  }

  InterpolateLinearArray() {}

  template <typename Iterator>
  void eval(T z, Iterator dest) const {
    const auto i = internal::interpolate_search(z, t_, false);
    if (i == t_.size() - 1) {
      std::copy_n(y_.begin() + i * n_, n_, dest);
    }
    const T r = (z - t_[i]) / (t_[i + 1] - t_[i]);
    auto y0 = y_.begin() + i * n_;
    auto y1 = y0 + n_;
    for (size_t j = 0; j < n_; ++j, ++y0, ++y1, ++dest) {
      *dest = *y0 + (*y1 - *y0) * r;
    }
  }
};

template <typename T>
class InterpolateSpline {
private:
  std::vector<T> t_;
  std::vector<T> y_;
  std::vector<T> k_;

public:
  InterpolateSpline(const std::vector<T> t, const std::vector<T>& y,
                    const char * name_t, const char * name_y) :
    t_(t), y_(y) {
    internal::validate_time(t_, y_.size(), name_t, name_y);
    const auto a = spline::calculate_a<T>(t);
    auto b = spline::calculate_b<T>(t, y);
    spline::solve_tridiagonal(a, b);
    k_ = b;
  }

  InterpolateSpline() {}

  T eval(T z) const {
    size_t i = internal::interpolate_search(z, t_, false);
    if (i == t_.size() - 1) {
      return y_[i];
    }
    const T t0 = t_[i], t1 = t_[i + 1], y0 = y_[i], y1 = y_[i + 1];
    const T k0 = k_[i], k1 = k_[i + 1];
    const T r = (z - t0) / (t1 - t0);
    const T a =  k0 * (t1 - t0) - (y1 - y0);
    const T b = -k1 * (t1 - t0) + (y1 - y0);
    return (1 - r) * y0 + r * y1 + r * (1 - r) * (a * (1 - r) + b * r);
  }
};

template <typename T, size_t rank>
class InterpolateSplineArray {
private:
  std::vector<T> t_;
  std::vector<T> y_;
  std::vector<T> k_;
  size_t n_;
public:
  InterpolateSplineArray(const std::vector<T> t, const std::vector<T>& y,
                           const array::dimensions<rank>& dim,
                           const char * name_t, const char * name_y) :
    t_(t), y_(y), k_(y_.size()), n_(dim.size) {
    const size_t len_y = y_.size() / n_;
    internal::validate_time(t, len_y, name_t, name_y);

    const auto a = spline::calculate_a<T>(t);
    std::vector<T> yi(t_.size());
    for (size_t i = 0; i < n_; ++i) {
      for (size_t j = 0; j < t_.size(); ++j) {
        yi[j] = y[i + j * n_];
      }
      auto b = spline::calculate_b<T>(t, yi);
      spline::solve_tridiagonal(a, b);
      for (size_t j = 0; j < t_.size(); ++j) {
        k_[i + j * n_] = b[j];
      }
    }
  }

  InterpolateSplineArray() {}

  template <typename Iterator>
  void eval(T z, Iterator dest) const {
    const auto i = internal::interpolate_search(z, t_, false);
    if (i == t_.size() - 1) {
      std::copy_n(y_.begin() + i * n_, n_, dest);
    }
    const auto t0 = t_[i];
    const auto t1 = t_[i + 1];
    const T r = (z - t0) / (t1 - t0);
    auto y0 = y_.begin() + i * n_;
    auto y1 = y0 + n_;
    auto k0 = k_.begin() + i * n_;
    auto k1 = k0 + n_;
    for (size_t j = 0; j < n_; ++j, ++y0, ++y1, ++k0, ++k1, ++dest) {
      const T a =  *k0 * (t1 - t0) - (*y1 - *y0);
      const T b = -*k1 * (t1 - t0) + (*y1 - *y0);
      *dest = (1 - r) * *y0 + r * *y1 + r * (1 - r) * (a * (1 - r) + b * r);
    }
  }
};

}
}
