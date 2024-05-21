#pragma once

#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace dust2 {
namespace details {

template <typename real_type>
real_type scale_log_weights(size_t n, typename std::vector<real_type>::iterator w) {
  if (n == 1) {
    return *w;
  }
  constexpr auto neg_infinity = -std::numeric_limits<real_type>::infinity();
  real_type max_w = neg_infinity;
  auto wi = w;
  for (size_t i = 0; i < n; ++i, ++wi) {
    if (std::isnan(*wi)) {
      *wi = neg_infinity;
    } else {
      max_w = std::max(max_w, *wi);
    }
  }

  if (max_w == neg_infinity) {
    return max_w;
  }

  real_type tot = 0.0;
  wi = w;
  for (size_t i = 0; i < n; ++i, ++wi) {
    *wi = std::exp(*wi - max_w);
    tot += *wi;
  }

  return std::log(tot / n) + max_w;
}


// This is a nasty bit of bookkeeping for the "systematic" resample,
// which is one of several options (we may allow configuration of this
// later, which will present its own challenges of course).
template <typename real_type>
void resample_weight(size_t n_particles,
                     typename std::vector<real_type>::const_iterator w,
                     const real_type u,
                     typename std::vector<size_t>::iterator idx) {
  const auto tot = std::accumulate(w, w + n_particles,
                                   static_cast<real_type>(0));
  const auto uu0 = tot * u / n_particles;
  const auto du = tot / n_particles;

  real_type ww = 0.0;
  size_t j = 0;
  for (size_t i = 0; i < n_particles; ++i) {
    const real_type uu = uu0 + i * du;
    // The second clause (i.e., j < n_particles) should never be hit but
    // prevents any invalid read if we have pathalogical 'u' that is
    // within floating point eps of 1 - solve this instead by passing
    // w as begin/end pair, something that only happens in single precision:
    //
    // https://github.com/reside-ic/reside-ic.github.io/pull/60/files
    // https://github.com/mrc-ide/dust/pull/238
    while (ww < uu && j < n_particles) {
      ww += *w;
      ++w;
      ++j;
    }
    *idx = j == 0 ? 0 : j - 1;
    ++idx;
  }
}


}
}
