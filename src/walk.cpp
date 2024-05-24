// Generated by dust2 (version 0.1.0) - do not edit
#include "../inst/examples/walk.cpp"

#include <cpp11.hpp>
#include <dust2/r/cpu.hpp>

[[cpp11::register]]
SEXP dust2_cpu_walk_alloc(cpp11::list r_pars,
                          cpp11::sexp r_time,
                          cpp11::sexp r_dt,
                          cpp11::sexp r_n_particles,
                          cpp11::sexp r_n_groups,
                          cpp11::sexp r_seed,
                          cpp11::sexp r_deterministic) {
  return dust2::r::dust2_cpu_alloc<walk>(r_pars, r_time, r_dt,
                                         r_n_particles, r_n_groups,
                                         r_seed, r_deterministic);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_run_steps(cpp11::sexp ptr, cpp11::sexp r_n_steps) {
  return dust2::r::dust2_cpu_run_steps<walk>(ptr, r_n_steps);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_cpu_run_to_time<walk>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_state(cpp11::sexp ptr, bool grouped) {
  return dust2::r::dust2_cpu_state<walk>(ptr, grouped);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_time(cpp11::sexp ptr) {
  return dust2::r::dust2_cpu_time<walk>(ptr);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_set_state_initial(cpp11::sexp ptr) {
  return dust2::r::dust2_cpu_set_state_initial<walk>(ptr);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_set_state(cpp11::sexp ptr, cpp11::sexp r_state,
                              bool grouped) {
  return dust2::r::dust2_cpu_set_state<walk>(ptr, r_state, grouped);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  return dust2::r::dust2_cpu_reorder<walk>(ptr, r_index);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_cpu_rng_state<walk>(ptr);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_cpu_set_time<walk>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_update_pars(cpp11::sexp ptr, cpp11::list pars,
                                bool grouped) {
  return dust2::r::dust2_cpu_update_pars<walk>(ptr, pars, grouped);
}

[[cpp11::register]]
SEXP dust2_cpu_walk_simulate(cpp11::sexp ptr, cpp11::sexp r_times,
                             cpp11::sexp r_index, bool grouped) {
  return dust2::r::dust2_cpu_simulate<walk>(ptr, r_times, r_index, grouped);
}
