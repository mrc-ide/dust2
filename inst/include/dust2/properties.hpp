#pragma once

#include <dust2/packing.hpp>
#include <dust2/continuous/delays.hpp>
#include <dust2/zero.hpp>
#include <vector>

namespace dust2 {

namespace internals {

template <class T, class = void>
struct test_has_build_internal: std::false_type {};
template <class T>
struct test_has_build_internal<T, std::void_t<decltype(T::build_internal)>>: std::true_type {};

template <class T, class = void>
struct test_has_packing_gradient: std::false_type {};
template <class T>
struct test_has_packing_gradient<T, std::void_t<decltype(T::packing_gradient)>>: std::true_type {};

template <class T, class = void>
struct test_has_zero_every: std::false_type {};
template <class T>
struct test_has_zero_every<T, std::void_t<decltype(T::zero_every)>>: std::true_type {};

template <class T, class = void>
struct test_has_update: std::false_type {};
template <class T>
struct test_has_update<T, std::void_t<decltype(T::update)>>: std::true_type {};

template <class T, class = void>
struct test_has_rhs: std::false_type {};
template <class T>
struct test_has_rhs<T, std::void_t<decltype(T::rhs)>>: std::true_type {};

template <class T, class = void>
struct test_has_output: std::false_type {};
template <class T>
struct test_has_output<T, std::void_t<decltype(T::output)>>: std::true_type {};

template <class T, class = void>
struct test_has_delays: std::false_type {};
template <class T>
struct test_has_delays<T, std::void_t<decltype(T::delays)>>: std::true_type {};

// These test that the signature of rhs and output consume the delays
// argument. Not especially lovely to read!
template <typename T>
using test_rhs_signature_uses_delays = std::is_invocable<decltype(T::rhs), typename T::real_type, typename T::real_type*, typename T::shared_state&, typename T::internal_state&, typename dust2::ode::delay_result_type<typename T::real_type>&, typename T::real_type*>;

template <typename T>
using test_output_signature_uses_delays = std::is_invocable<decltype(T::output), typename T::real_type, typename T::real_type*, typename T::shared_state&, typename T::internal_state&, typename dust2::ode::delay_result_type<typename T::real_type>&>;

// Because not all functions have output, and because I can't work out
// how to get template checks to short circuit, I have done this with
// 'if constexpr' which *does* short circuit, getting us a compile
// time result that checks for all of:
// * existance of delays
// * existance of output
// * an output function that consumes the delay
template <typename T>
constexpr bool test_output_uses_delays() {
  if constexpr (test_has_delays<T>::value && test_has_output<T>::value) {
    return test_output_signature_uses_delays<T>::value;
  } else {
    return false;
  }
}

template <typename T>
constexpr bool test_rhs_uses_delays() {
  if constexpr (test_has_delays<T>::value) {
    return test_rhs_signature_uses_delays<T>::value;
  } else {
    return false;
  }
}

}

template<typename T>
struct properties {
  using has_build_internal = internals::test_has_build_internal<T>;
  using has_packing_gradient = internals::test_has_packing_gradient<T>;
  using has_zero_every = internals::test_has_zero_every<T>;
  using is_mixed_time = typename std::conditional<internals::test_has_rhs<T>::value && internals::test_has_update<T>::value, std::true_type, std::false_type>::type;
  using has_output = typename std::conditional<internals::test_has_rhs<T>::value && internals::test_has_output<T>::value, std::true_type, std::false_type>::type;
  using has_delays = internals::test_has_delays<T>;
  // Because of the above these are now actual numbers rather than types; we may make this change everywhere...
  static constexpr bool rhs_uses_delays = internals::test_rhs_uses_delays<T>();
  static constexpr bool output_uses_delays = internals::test_output_uses_delays<T>();
};

// wrappers around some uses of member functions that may or may not
// exist, centralising most of the weird into this file:
template <typename T, typename std::enable_if<properties<T>::has_build_internal::value, T>::type* = nullptr>
typename T::internal_state do_build_internal(const typename T::shared_state &shared) {
  return T::build_internal(shared);
}

template <typename T, typename std::enable_if<!properties<T>::has_build_internal::value, T>::type* = nullptr>
typename T::internal_state do_build_internal(const typename T::shared_state &shared) {
  return typename T::internal_state{};
}

template <typename T, typename std::enable_if<properties<T>::has_packing_gradient::value, T>::type* = nullptr>
dust2::packing do_packing_gradient(const typename T::shared_state& shared) {
  return T::packing_gradient(shared);
}

template <typename T, typename std::enable_if<!properties<T>::has_packing_gradient::value, T>::type* = nullptr>
dust2::packing do_packing_gradient(const typename T::shared_state &shared) {
  return dust2::packing{};
}


template <typename T, typename std::enable_if<properties<T>::has_output::value, T>::type* = nullptr>
size_t do_n_state_output(const dust2::packing& packing) {
  return std::accumulate(packing.len().end() - T::size_output(),
                         packing.len().end(),
                         0);
}

template <typename T, typename std::enable_if<!properties<T>::has_output::value, T>::type* = nullptr>
size_t do_n_state_output(const dust2::packing& packing) {
  return 0;
}

template <typename T, typename std::enable_if<properties<T>::has_zero_every::value, T>::type* = nullptr>
auto zero_every_vec(const std::vector<typename T::shared_state>& shared) {
  using real_type = typename T::real_type;
  std::vector<zero_every_type<real_type>> ret;
  ret.reserve(shared.size());
  for (const auto& el : shared) {
    ret.push_back(T::zero_every(el));
  }
  return ret;
}

template <typename T, typename std::enable_if<!properties<T>::has_zero_every::value, T>::type* = nullptr>
auto zero_every_vec(const std::vector<typename T::shared_state>& shared) {
  using real_type = typename T::real_type;
  return std::vector<zero_every_type<real_type>>(shared.size(), dust2::zero_every_type<real_type>());
}

template <typename T, typename std::enable_if<properties<T>::has_delays::value, T>::type* = nullptr>
auto do_delays(const std::vector<typename T::shared_state>& shared) {
  using real_type = typename T::real_type;
  std::vector<ode::delays<real_type>> ret;
  ret.reserve(shared.size());
  for (size_t i = 0; i < shared.size(); ++i) {
    ret.push_back(T::delays(shared[i]));
  }
  return ret;
}

template <typename T, typename std::enable_if<!properties<T>::has_delays::value, T>::type* = nullptr>
auto do_delays(const std::vector<typename T::shared_state>& shared) {
  using real_type = typename T::real_type;
  return std::vector<ode::delays<real_type>>(shared.size(), dust2::ode::delays<real_type>{{}});
}

}
