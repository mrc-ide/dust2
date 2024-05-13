 #pragma once

#include <dust2/r/helpers.hpp>
#include <dust2/cpu.hpp>

namespace dust2 {
namespace r {

// There's a big question here of if we take a model pointer or do the
// initialisation ourself.  I think that we really need to do the
// latter because otherwise the rng is quite hard to think about
// because the user might end up holding two copies?  Or we can offer
// both approaches (though we've never once needed to initialise the
// model first).  The user will mostly interact with this as
//
// filter <- make_filter_sir(pars, ...)
// filter$run(pars) # update pars, verify size has not changed.
//
// The first bit of init
template <typename T>
struct unfilter_state {
  using real_type = typename T::real_type;
  using data_type = typename T::data_type;
  dust_cpu<T> model;
  real_type time_start;
  std::vector<real_type> time;
  std::vector<size_t> step;
  std::vector<std::vector<data_type>> data;
  size_t n_groups;
};

// TODO: this name must be changed!
template <typename T>
cpp11::sexp dust2_cpu_unfilter_alloc(cpp11::list r_pars,
                                     cpp11::sexp r_time_start,
                                     cpp11::sexp r_time,
                                     cpp11::sexp r_dt,
                                     cpp11::list r_data,
                                     cpp11::sexp r_n_groups) {
  using real_type = typename T::real_type;
  using rng_state_type = typename T::rng_state_type;

  const auto n_groups = to_size(r_n_groups, "n_groups");
  const auto time_start = check_time(r_time_start, "time_start");
  const auto time = check_time_sequence<real_type>(time_start, r_time, "time");
  const auto dt = check_dt(r_dt);
  const auto shared = build_shared<T>(r_pars, n_groups);
  const auto internal = build_internal<T>(shared);
  const auto data = check_data<T>(r_data, time.size(), n_groups, "data");

  // This probably gets reused a bit, too, easily pulled into a
  // helper; it will be used in a simulate() method too
  std::vector<size_t> step;
  for (size_t i = 0; i < time.size(); i++) {
    const auto t0 = i == 0 ? time_start : time[i - 1];
    const auto t1 = time[i];
    step.push_back(static_cast<size_t>(std::round((t1 - t0) / dt)));
  }

  // It's possible that we don't want to always really be
  // deterministic here?  Though nooone can think of a case where
  // that's actually the behaviour wanted.  For now let's go fully
  // deterministic.
  auto seed = mcstate::random::r::as_rng_seed<rng_state_type>(R_NilValue);
  const auto deterministic = true;
  const size_t n_particles = 1;

  // Then allocate the model; this pulls together almost all the data
  // we need.  At this point we could have constructed the model out
  // of one that exists already on the R side, but I think that's
  // going to feel weirder overall.
  const auto model = dust_cpu<T>(shared, internal, time_start, dt, n_particles,
                                 seed, deterministic);

  auto obj = new unfilter_state<T>{model, time_start, time, step, data,
                                   shared.size()};
  cpp11::external_pointer<unfilter_state<T>> ptr(obj, true, false);

  // NOTE: not (yet) returning length here, but we could
  const bool grouped = n_groups > 0;
  cpp11::sexp r_n_state = cpp11::as_sexp(obj->model.n_state());
  cpp11::sexp r_group_names = R_NilValue;
  //cpp11::sexp r_group_names = grouped ? r_pars.attr("names") : cpp11::as_sexp(R_NilValue);
  cpp11::sexp r_grouped = cpp11::as_sexp(grouped);

  return cpp11::writable::list{ptr, r_n_state, r_grouped, r_group_names};
}

template <typename T>
cpp11::sexp dust2_cpu_unfilter_run(cpp11::sexp ptr, cpp11::sexp r_pars) {
  using real_type = typename T::real_type;
  auto *obj =
    cpp11::as_cpp<cpp11::external_pointer<unfilter_state<T>>>(ptr).get();

  if (r_pars != R_NilValue) {
    // obj->model->update_pars - pull out the parameter handling code here.
  }

  // There's a question of how much we move this into the C++ part
  // (rather than this interface layer) which we will want to do to
  // make this more extendable.
  const auto time_start = obj->time_start;
  const auto& step = obj->step;
  const auto n_groups = obj->n_groups;
  const auto n_times = step.size();

  // Do we want these stored into the object? Seems a low cost really?
  std::vector<real_type> ll_step(n_groups, 0);
  std::vector<real_type> ll(n_groups, 0);

  // obj->model.update_pars(pars);
  obj->model.set_time(time_start);
  obj->model.set_state_initial();

  for (size_t i = 0; i < n_times; ++i) {
    obj->model.run_steps(step[i]); // just compute this at point of use?
    obj->model.compare_data(obj->data[i], ll_step.begin());
    for (size_t j = 0; j < n_groups; ++j) {
      ll[j] += ll_step[j];
    }
  }

  cpp11::writable::doubles r_ll(ll);

  return r_ll;
}


}
}
