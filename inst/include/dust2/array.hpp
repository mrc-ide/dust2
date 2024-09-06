#pragma once

#include <cstddef>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace dust2 {
namespace array {

template <size_t rank_>
struct dimensions {
  size_t rank;
  std::array<size_t, rank_> dim;
  size_t size;
  std::array<size_t, rank_> mult;
  // One disadvantage here is we can't optimise away the case where
  // arrays have literal size, we might want to add that in as an
  // optimisation later?
  dimensions(std::initializer_list<size_t> dim_) :
    rank(rank_),
    dim(dim_),
    size(std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<size_t>())) {
    // TODO: I think we can fail to compile here if rank does not
    // match dim_.size()?
    //
    // TODO: I think we can drop 2 of these, so only bother doing
    // this in the case where we have 3 dimensions or more and store 2
    // fewer values (with slightly more complex calculations).
    for (size_t i = 0; i < rank; ++i) {
      mult[i] = i == 0 ? 1 : dim[i - 1] * mult[i - 1];
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

size_t index(std::initializer_list<size_t> at,
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
