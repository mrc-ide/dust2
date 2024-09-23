#pragma once

#include <array>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <monty/random/math.hpp>

namespace dust2 {
namespace array {

namespace {

constexpr auto min2 = [](auto a, auto b) { return a < b ? a : b; };
constexpr auto max2 = [](auto a, auto b) { return a > b ? a : b; };

}

template <size_t rank_>
struct dimensions {
  size_t rank;
  size_t size;
  std::array<size_t, rank_> dim;
  std::array<size_t, rank_> mult;

  template <typename Iterator>
  dimensions(Iterator iter) : rank(rank_), size(1) {
    for (size_t i = 0; i < rank; ++i, ++iter) {
      dim[i] = *iter;
      mult[i] = size;
      size *= dim[i];
    }
  }

  dimensions(std::initializer_list<size_t> dim_) : dimensions(dim_.begin()) {
  }
};

inline size_t index1(size_t at, size_t size, const char *name,
                     const char *context) {
  if (at >= size) {
    std::stringstream msg;
    msg << "Tried to access element " << at + 1 << " of vector " << name <<
      ", which is out of bounds (length is " << size << ").";
    if (context != nullptr && *context != 0) {
      msg << " " << context;
    }
    throw std::runtime_error(msg.str());
  }
  return at;
}

inline size_t index2(size_t at_r, size_t at_c, size_t size_r, size_t size_c,
                     const char *name, const char *context) {
  if (at_r >= size_r || at_c >= size_c) {
    std::stringstream msg;
    msg << "Tried to access element (" << at_r + 1 << ", " << at_c + 1 <<
      ") of matrix " << name <<
      ", which is out of bounds (" << size_r << "rows and " << size_c <<
      "columns).";
    if (context != nullptr && *context != 0) {
      msg << " " << context;
    }
    throw std::runtime_error(msg.str());
  }
  return at_r + at_c * size_r;
}

template <size_t rank>
size_t index(std::initializer_list<size_t> at,
             const dimensions<rank> dim,
             const char *name,
             const char *context) {
  size_t ret = 0;
  auto it = at.begin();
  for (size_t i = 0; i < rank; ++i, ++it) {
    const auto v = *it;
    if (v >= dim.dim[i]) {
      std::stringstream msg;
      msg << "Tried to access element " << v + 1 << " from dimension " <<
        i + 1 << " of array " << name <<
        ", which is out of bounds (extent of this dimension is " <<
        dim.dim[i] << ").";
      if (context != nullptr && *context != 0) {
        msg << " " << context;
      }
      throw std::runtime_error(msg.str());
    }
    ret += dim.mult[i] * v;
  }
  return ret;
}

namespace {
using idx = std::pair<size_t, size_t>;
}

// Sum over vector (1d)
template <typename T, typename Container, typename BinaryOp>
T reduce(Container x, const dimensions<1>& dim, T init, BinaryOp op,
         const idx& i) {
  if (i.first > i.second) {
    return init;
  }
  return std::accumulate(x + i.first, x + i.second + 1, init, op);
}

// These are written out by hand for now, and then we can explore less
// tedious ways, and for more dimensions.  Previously we've generated
// code for this rather than trying to do anything too clever with C++
// template metaprogramming.

// 2 dimensions:
template <typename T, typename Container, typename BinaryOp>
T reduce(Container x, const dimensions<2>& dim, T init, BinaryOp op,
         const idx& i, const idx& j) {
  T tot = init;
  for (size_t jj = j.first; jj <= j.second; ++jj) {
    for (size_t ii = i.first; ii <= i.second; ++ii) {
      tot = op(tot, x[ii + jj * dim.mult[1]]);
    }
  }
  return tot;
}

// 3 dimensions:
template <typename T, typename Container, typename BinaryOp>
T reduce(Container x, const dimensions<3>& dim, T init, BinaryOp op,
         const idx& i, const idx& j, const idx& k) {
  T tot = init;
  for (size_t kk = k.first; kk <= k.second; ++kk) {
    for (size_t jj = j.first; jj <= j.second; ++jj) {
      for (size_t ii = i.first; ii <= i.second; ++ii) {
        tot = op(tot, x[ii + jj * dim.mult[1] + kk * dim.mult[2]]);
      }
    }
  }
  return tot;
}

// Special case, sum over everything:
template <typename T, typename Container, typename BinaryOp, size_t rank>
T reduce(Container x, const dimensions<rank>& dim, T init, BinaryOp op) {
  return reduce(x, x + dim.size, init, op);
}

// Then implement sums:
template <typename T, typename Container, size_t rank>
T sum(Container x, const dimensions<rank>& dim) {
  return reduce(x, dim, static_cast<T>(0), std::plus<T>());
}

template <typename T, typename Container>
T sum(Container x, const dimensions<1>& dim, const idx& i) {
  return reduce(x, dim, static_cast<T>(0), std::plus<T>(), i);
}

template <typename T, typename Container>
T sum(Container x, const dimensions<2>& dim, const idx& i, const idx& j) {
  return reduce(x, dim, static_cast<T>(0), std::plus<T>(), i, j);
}

template <typename T, typename Container>
T sum(Container x, const dimensions<3>& dim, const idx& i, const idx& j, const idx& k) {
  return reduce(x, dim, static_cast<T>(0), std::plus<T>(), i, j, k);
}

// And products:
template <typename T, typename Container, size_t rank>
T prod(Container x, const dimensions<rank>& dim) {
  return reduce(x, dim, static_cast<T>(1), std::multiplies<T>());
}

template <typename T, typename Container>
T prod(Container x, const dimensions<1>& dim, const idx& i) {
  return reduce(x, dim, static_cast<T>(1), std::multiplies<T>(), i);
}

template <typename T, typename Container>
T prod(Container x, const dimensions<2>& dim, const idx& i, const idx& j) {
  return reduce(x, dim, static_cast<T>(1), std::multiplies<T>(), i, j);
}

template <typename T, typename Container>
T prod(Container x, const dimensions<3>& dim, const idx& i, const idx& j, const idx& k) {
  return reduce(x, dim, static_cast<T>(1), std::multiplies<T>(), i, j, k);
}

// And min:
template <typename T, typename Container, size_t rank>
T min(Container x, const dimensions<rank>& dim) {
  // TODO: I don't see why I can't use
  // > return reduce(x, dim, std::numeric_limits<T>::infinity(), min2);
  // but I get a template substitution error on compilation if I do so.
  return std::accumulate(x, x + dim.size, std::numeric_limits<T>::infinity(), min2);
}

template <typename T, typename Container>
T min(Container x, const dimensions<1>& dim, const idx& i) {
  return reduce(x, dim, std::numeric_limits<T>::infinity(), min2, i);
}

template <typename T, typename Container>
T min(Container x, const dimensions<2>& dim, const idx& i, const idx& j) {
  return reduce(x, dim, std::numeric_limits<T>::infinity(), min2, i, j);
}

template <typename T, typename Container>
T min(Container x, const dimensions<3>& dim, const idx& i, const idx& j, const idx& k) {
  return reduce(x, dim, std::numeric_limits<T>::infinity(), min2, i, j, k);
}

// And max
template <typename T, typename Container, size_t rank>
T max(Container x, const dimensions<rank>& dim) {
  return std::accumulate(x, x + dim.size, std::numeric_limits<T>::infinity(), max2);
}

template <typename T, typename Container>
T max(Container x, const dimensions<1>& dim, const idx& i) {
  return reduce(x, dim, -std::numeric_limits<T>::infinity(), max2, i);
}

template <typename T, typename Container>
T max(Container x, const dimensions<2>& dim, const idx& i, const idx& j) {
  return reduce(x, dim, -std::numeric_limits<T>::infinity(), max2, i, j);
}

template <typename T, typename Container>
T max(Container x, const dimensions<3>& dim, const idx& i, const idx& j, const idx& k) {
  return reduce(x, dim, -std::numeric_limits<T>::infinity(), max2, i, j, k);
}

}
}
