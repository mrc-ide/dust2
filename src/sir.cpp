// Generated by dust2 (version 0.1.0) - do not edit
#include <cpp11.hpp>
#include <dust2/r/cpu.hpp>
#include <dust2/r/filter.hpp>

// first declarations all go here, with their decorators, once we get
// this bit sorted.

#include "../inst/examples/sir.cpp"

[[cpp11::register]]
SEXP dust2_cpu_sir_alloc(cpp11::list r_pars,
                         cpp11::sexp r_time,
                         cpp11::sexp r_dt,
                         cpp11::sexp r_n_particles,
                         cpp11::sexp r_n_groups,
                         cpp11::sexp r_seed,
                         cpp11::sexp r_deterministic) {
  return dust2::r::dust2_cpu_alloc<sir>(r_pars, r_time, r_dt,
                                        r_n_particles, r_n_groups,
                                        r_seed, r_deterministic);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_run_steps(cpp11::sexp ptr, cpp11::sexp r_n_steps) {
  return dust2::r::dust2_cpu_run_steps<sir>(ptr, r_n_steps);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_state(cpp11::sexp ptr, bool grouped) {
  return dust2::r::dust2_cpu_state<sir>(ptr, grouped);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_time(cpp11::sexp ptr) {
  return dust2::r::dust2_cpu_time<sir>(ptr);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_update_pars(cpp11::sexp ptr, cpp11::list pars,
                               bool grouped) {
  return dust2::r::dust2_cpu_update_pars<sir>(ptr, pars, grouped);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_set_state_initial(cpp11::sexp ptr) {
  return dust2::r::dust2_cpu_set_state_initial<sir>(ptr);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_set_state(cpp11::sexp ptr, cpp11::sexp r_state) {
  return dust2::r::dust2_cpu_set_state<sir>(ptr, r_state);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_cpu_rng_state<sir>(ptr);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_compare_data(cpp11::sexp ptr,
                                cpp11::sexp r_data,
                                bool grouped) {
  return dust2::r::dust2_cpu_compare_data<sir>(ptr, r_data, grouped);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_unfilter_alloc(cpp11::list r_pars,
                                  cpp11::sexp r_time_start,
                                  cpp11::sexp r_time,
                                  cpp11::sexp r_dt,
                                  cpp11::list r_data,
                                  cpp11::sexp r_n_groups) {
  return dust2::r::dust2_cpu_unfilter_alloc<sir>(r_pars, r_time_start, r_time,
                                                 r_dt, r_data, r_n_groups);
}

[[cpp11::register]]
SEXP dust2_cpu_sir_unfilter_run(cpp11::sexp ptr, cpp11::sexp r_pars,
                                bool grouped) {
  return dust2::r::dust2_cpu_unfilter_run<sir>(ptr, r_pars, grouped);
}
