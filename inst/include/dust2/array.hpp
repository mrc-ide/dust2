#pragma once

#include <array>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace dust2 {
namespace array {

template <size_t rank_>
struct dimensions {
  size_t rank;
  size_t size;
  std::array<size_t, rank_> dim;
  std::array<size_t, rank_> mult;
  // One disadvantage here is we can't optimise away the case where
  // arrays have literal size, we might want to add that in as an
  // optimisation later?
  dimensions(std::initializer_list<size_t> dim_) :
    rank(rank_),
    size(1) {
    auto iter = dim_.begin();
    for (size_t i = 0; i < rank; ++i, ++iter) {
      dim[i] = *iter;
      mult[i] = size;
      size *= dim[i];
    }
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

// One dimension:
template <typename T, typename Iterator>
T sum(const Iterator& x, const idx& i) {
  return std::accumulate(x + i.first, x + i.second, static_cast<T>(0));
}

// Complete sum, a special case of 1d but applicable to any rank:
template <typename T, typename Container, size_t rank>
T sum(const Container& x, const dimensions<rank>& dim) {
  return sum<T>(x, {static_cast<size_t>(0), dim.size});
}

// These are written out by hand for now, and then we an explore less
// tedious ways, and for more dimensions.  I think we can do this
// really quite generally with iterators?
// two dimensions:
template <typename T, typename Container>
T sum(Container x, const idx& i, const idx& j, const dimensions<2>& dim) {
  T tot = 0;
  for (int jj = j.first; jj < j.second; ++jj) {
    for (int ii = i.first; ii < i.second; ++ii) {
      tot += x[ii + jj * dim.mult[1]];
    }
  }
  return tot;
}

/*
// three dimensions:
template <typename T, typename Container>
T sum(const Container& x, const idx& i, const idx& j, const idx& k, const dimensions<3>& dim) {
  T tot = 0;
  for (int kk = k.first; kk < k.second; ++kk) {
    for (int jj = j.first; jj < j.second; ++jj) {
      for (int ii = i.first; ii < i.second; ++ii) {
        tot += x[ii + jj * dim.mult[1] + kk * dim.mult[2]];
      }
    }
  }
  return tot;
 }
*/

}
}
