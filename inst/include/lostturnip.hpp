#pragma once
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lostturnip {

// Declaring these here, rather than within the find_result, as
// otherwise we get a compiler warning about using experimental cuda
// features. It will be equivalent though, but does require C++14.
namespace {
template <typename real_type>
constexpr real_type na = std::numeric_limits<real_type>::quiet_NaN();

template <typename real_type>
constexpr real_type eps = std::numeric_limits<real_type>::epsilon();
}

template <typename real_type>
struct result {
  real_type x;
  real_type fx;
  int iterations;
  bool converged;
};

// From zeroin.c, in brent.shar
template <typename real_type, typename F>
#ifdef __NVCC__
__host__ __device__
#endif
result<real_type> find_result(F f, real_type a, real_type b,
                              real_type tol, int max_iterations) {
  real_type fa = f(a);
  real_type fb = f(b);
  int iterations = 0;
  bool converged = false;

  if (fa == 0) {
    b = a;
    fb = fa;
    converged = true;
  } else if (fb == 0) {
    converged = true;
  } else if (fa * fb > 0) {
    // Same sign; can't find root with this:
    b = na<real_type>;
    fb = na<real_type>;
    converged = false;
  } else {
    real_type c = a;
    real_type fc = fa;   // c = a, f(c) = f(a)

    for (; iterations < max_iterations; ++iterations) { // Main iteration loop
      // Distance from the last but one to the last approximation
      const real_type prev_step = b - a;

      // Interpolation step is calculated in the form p/q; division
      // operations is dlayed until the last moment
      real_type p;
      real_type q;

      if (std::abs(fc) < std::abs(fb)) {
        // Swap data for b to be the best approximation
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
      }

      // Actual tolerance
      const real_type tol_act = 2 * eps<real_type> * std::abs(b) + tol / 2;
      // Step at this iteration
      real_type new_step = (c - b) / 2;

      if (std::abs(new_step) <= tol_act || fb == 0) {
        // Acceptable approximation is found
        converged = true;
        break;
      }

      // increase readability below, avoids many repeated static casts
      const real_type one = 1;

      // Decide if the interpolation can be tried
      //
      // If prev_step was large enough and was in true direction, then
      // interpolation can be tried
      if (std::abs(prev_step) >= tol_act && std::abs(fa) > std::abs(fb)) {
        // interpolation
        const real_type cb = c - b;
        if (a == c) {
          // If we have only two distinct points linear interpolation
          // can only be applied
          const real_type t1 = fb / fa;
          p = cb * t1;
          q = one - t1;
        } else {
          // Quadric inverse interpolation
          q = fa / fc;
          const real_type t1 = fb / fc;
          const real_type t2 = fb / fa;
          p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - one));
          q = (q - one) * (t1 - one) * (t2 - one);
        }
        if (p > 0) {
          // p was calculated with the opposite sign; make p positive
          // and assign possible minus to q
          q = -q;
        } else {
          p = -p;
        }

        // If b + p / q falls in [b, c] and isn't too large it is
        // accepted
        //
        // If p / q is too large then the bissection procedure can
        // reduce [b,c] range to more extent
        if (p < (static_cast<real_type>(0.75) * cb * q - std::abs(tol_act * q) / 2) &&
            p < std::abs(prev_step * q / 2)) {
          new_step = p / q;
        }
      }

      // Adjust the step to be not less than tolerance
      if (std::abs(new_step) < tol_act) {
        new_step = std::copysign(tol_act, new_step);
      }

      // Save the previous approximation
      a = b;
      fa = fb;
      // Do step to a new approximation
      b += new_step;
      fb = f(b);
      if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
        // Adjust c for it to have a sign opposite to that of b
        c = a;  fc = fa;
      }
    }
  }

#ifdef __CUDA_ARCH__
  __syncwarp();
#endif
  return result<real_type>{b, fb, iterations, converged};
}

template <typename real_type, typename F>
#ifdef __NVCC__
__host__ __device__
#endif
real_type find(F f, real_type a, real_type b,
               real_type tol, int max_iterations) {
  const auto result = find_result(f, a, b, tol, max_iterations);
  if (!result.converged) {
#ifdef __CUDA_ARCH__
    printf("some error\n");
    __trap();
#else
    throw std::runtime_error("some error");
#endif
  }
  return result.x;
}

}
