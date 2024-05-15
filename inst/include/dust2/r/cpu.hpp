#pragma once

#include <dust2/r/helpers.hpp>
#include <mcstate/r/random.hpp>

#include <dust2/cpu.hpp>

namespace dust2 {
namespace r {

template <typename T>
SEXP dust2_cpu_alloc(cpp11::list r_pars,
                     cpp11::sexp r_time,
                     cpp11::sexp r_dt,
                     cpp11::sexp r_n_particles,
                     cpp11::sexp r_n_groups,
                     cpp11::sexp r_seed,
                     cpp11::sexp r_deterministic) {
  using rng_state_type = typename T::rng_state_type;

  const auto time = check_time(r_time, "time");
  const auto dt = check_dt(r_dt);

  auto n_particles = to_size(r_n_particles, "n_particles");
  auto n_groups = to_size(r_n_groups, "n_groups");

  const auto shared = build_shared<T>(r_pars, n_groups);
  // Later, we need one of these per thread
  const auto internal = build_internal<T>(shared);

  auto seed = mcstate::random::r::as_rng_seed<rng_state_type>(r_seed);
  auto deterministic = to_bool(r_deterministic, "deterministic");

  auto obj = new dust_cpu<T>(shared, internal, time, dt, n_particles,
                             seed, deterministic);
  cpp11::external_pointer<dust_cpu<T>> ptr(obj, true, false);

  // Later, we'll export a bit more back from the model (in particular
  // models need to provide information about how they organise
  // variables, ode models report computed control, etc.
  const auto grouped = n_groups > 0;
  cpp11::sexp r_group_names = R_NilValue;
  if (grouped) {
    r_group_names = r_pars.attr("names");
  }
  cpp11::sexp r_n_state = cpp11::as_sexp(obj->n_state());
  cpp11::sexp r_grouped = cpp11::as_sexp(grouped);

  return cpp11::writable::list{ptr, r_n_state, r_grouped, r_group_names};
}

template <typename T>
SEXP dust2_cpu_run_steps(cpp11::sexp ptr, cpp11::sexp r_n_steps) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  auto n_steps = to_size(r_n_steps, "n_steps");
  obj->run_steps(n_steps);
  return R_NilValue;
}

template <typename T>
SEXP dust2_cpu_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  const auto time = check_time(r_time, "time");
  const auto curr = obj->time();
  if (time < curr) {
    cpp11::stop("Can't run to time %f, model already at time %f",
                time, curr);
  }
  obj->run_to_time(time);
  return R_NilValue;
}

template <typename T>
SEXP dust2_cpu_state(cpp11::sexp ptr, bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  cpp11::sexp ret = R_NilValue;
  const auto iter = obj->state().begin();
  if (grouped) {
    ret = export_array_n(iter,
                         {obj->n_state(), obj->n_particles(), obj->n_groups()});
  } else {
    ret = export_array_n(iter,
                         {obj->n_state(), obj->n_particles() * obj->n_groups()});
  }
  return ret;
}

template <typename T>
SEXP dust2_cpu_time(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  return cpp11::as_sexp(obj->time());
}

// If this is a grouped model then we will return a matrix with
// dimensions )(state x particle x group) and if we are ungrouped
// (state x particle); the difference is really only apparent in the
// case where we have a single group to drop.  In dust1 we had the
// option for (state x group) for simulations too.
//
// For now perhaps lets just ignore this detail and always return the
// 3d version as the code to do grouped/ungrouped switching is deep
// within one of the many in-progress branches and I can't find it
// yet.
//
// In the case where we have time, that goes in the last position
template <typename T>
SEXP dust2_cpu_set_state_initial(cpp11::sexp ptr) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  obj->set_state_initial();
  return R_NilValue;
}

template <typename T>
SEXP dust2_cpu_set_state(cpp11::sexp ptr, cpp11::sexp r_state, bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  // Suppose that we have a n_state x n_particles x n_groups grouped
  // system, we then require that we have a state array with rank 3;
  // for an ungrouped system this will be rank 2 array.
  auto dim = cpp11::as_cpp<cpp11::integers>(r_state.attr("dim"));
  const auto rank = dim.size();
  const auto rank_expected = grouped ? 3 : 2;
  if (rank != rank_expected) {
    cpp11::stop("Expected 'state' to be a %dd array", rank_expected);
  }
  const int n_state = obj->n_state();
  const int n_particles =
    grouped ? obj->n_particles() : obj->n_particles() * obj->n_groups();
  const int n_groups = grouped ? obj->n_groups() : 1;
  if (dim[0] != n_state) {
    cpp11::stop("Expected the first dimension of 'state' to have size %d",
                n_state);
  }
  const auto recycle_particle = n_particles > 1 && dim[1] == 1;
  if (dim[1] != n_particles && dim[1] != 1) {
    cpp11::stop("Expected the second dimension of 'state' to have size %d or 1",
                n_particles);
  }
  const auto recycle_group = !grouped || (n_groups > 1 && dim[2] == 1);
  if (grouped && dim[2] != n_groups && dim[2] != 1) {
    cpp11::stop("Expected the third dimension of 'state' to have size %d or 1",
                n_groups);
  }
  obj->set_state(REAL(r_state), recycle_particle, recycle_group);
  return R_NilValue;
}

// Not really intended to be called by users, this just helps us test
// bookkeeping really.
template <typename T>
SEXP dust2_cpu_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  const auto n_particles = obj->n_particles();
  const auto n_groups = obj->n_groups();
  const int len = n_particles * n_groups;
  // We really should expect a matrix perhaps, but this test will
  // catch that confusingly at least.  Users won't actually call this.
  if (r_index.size() != len) {
    cpp11::stop("Expected an index of length %d", len);
  }
  std::vector<size_t> index;
  index.reserve(len);
  for (auto i : r_index) {
    if (i < 1 || i > len) {
      cpp11::stop("Expected 'index' values to lie in [1, %d]", n_particles);
    }
    index.push_back(i - 1);
  }
  obj->reorder(index.begin());
  return R_NilValue;
}

template <typename T>
SEXP dust2_cpu_rng_state(cpp11::sexp ptr) {
  using rng_state_type = typename T::rng_state_type;
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();

  const auto state = obj->rng_state();
  const auto len = sizeof(typename rng_state_type::int_type) * state.size();
  cpp11::writable::raws ret(len);
  std::memcpy(RAW(ret), state.data(), len);
  return ret;
}

template <typename T>
SEXP dust2_cpu_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  const auto time = check_time(r_time, "time");
  obj->set_time(time);
  return R_NilValue;
}

template <typename T>
SEXP dust2_cpu_update_pars(cpp11::sexp ptr, cpp11::list r_pars,
                           bool grouped) {
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  update_pars(*obj, r_pars, grouped);
  return R_NilValue;
}

// This one exists to help push around the comparison part of things;
// it's not expected to be called often by users.
template <typename T>
SEXP dust2_cpu_compare_data(cpp11::sexp ptr,
                            cpp11::sexp r_data,
                            bool grouped) {
  using data_type = typename T::data_type;
  auto *obj = cpp11::as_cpp<cpp11::external_pointer<dust_cpu<T>>>(ptr).get();
  const auto n_groups = obj->n_groups();
  std::vector<data_type> data;
  if (grouped) {
    auto r_data_list = cpp11::as_cpp<cpp11::list>(r_data);
    check_length(r_data_list, n_groups, "data");
    for (size_t i = 0; i < n_groups; ++i) {
      data.push_back(T::build_data(r_data_list[i]));
    }
  } else {
    if (n_groups > 1) {
      cpp11::stop("Can't compare with grouped = FALSE with more than one group");
    }
    data.push_back(T::build_data(r_data));
  }

  cpp11::writable::doubles ret(obj->n_particles() * obj->n_groups());
  obj->compare_data(data.begin(), REAL(ret));
  if (grouped) {
    set_array_dims(ret, {obj->n_particles(), obj->n_groups()});
  }
  return ret;
}

}
}
