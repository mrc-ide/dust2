#pragma once

#include <cstddef>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace dust2 {
namespace array {

struct array_dimensions {
  std::vector<size_t> dim;
  size_t rank;
  size_t size;
  std::vector<size_t> mult;
  // One disadvantage here is we can't optimise away the case where
  // arrays have literal size, we might want to add that in as an
  // optimisation later?
  array_dimensions(std::initializer_list<size_t> dim_) :
    dim(dim_),
    rank(dim.size()),
    size(std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<size_t>())) {
    mult.reserve(dim.size());
    for (size_t i = 0; i < rank; ++i) {
      mult.push_back(i == 0 ? 1 : dim[i - 1] * mult[i - 1]);
    }
  }
};

inline size_t position1(size_t at, size_t size, const char *name,
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

inline size_t position2(size_t at_r, size_t at_c, size_t size_r, size_t size_c,
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


size_t position(std::initializer_list<size_t> at,
                const std::vector<size_t>& dim,
                const std::vector<size_t>& mult,
                const char *name,
                const char *context) {
  const auto rank = dim.size();
  size_t ret = 0;
  auto it = at.begin();
  for (size_t i = 0; i < rank; ++i, ++it) {
    const auto v = *it;
    if (v >= dim[i]) {
      std::stringstream msg;
      msg << "Tried to access element " << v + 1 << " from dimension " <<
        i + 1 << " of array " << name <<
        ", which is out of bounds (extent of this dimension is " <<
        dim[i] << ").";
      if (context != nullptr && *context != 0) {
        msg << " " << context;
      }
      throw std::runtime_error(msg.str());
    }
    ret += mult[i] * v;
  }
  return ret;
}

}
}
