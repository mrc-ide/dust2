#pragma once

#include <dust2/packing.hpp>
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

}

template<typename T>
struct properties {
  using has_build_internal = internals::test_has_build_internal<T>;
  using has_packing_gradient = internals::test_has_packing_gradient<T>;
  using has_zero_every = internals::test_has_zero_every<T>;
  using is_mixed_time = typename std::conditional<internals::test_has_rhs<T>::value && internals::test_has_update<T>::value, std::true_type, std::false_type>::type;
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

}
