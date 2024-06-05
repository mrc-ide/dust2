[[cpp11::register]]
SEXP dust2_cpu_{{name}}_alloc(cpp11::list r_pars, cpp11::sexp r_time, cpp11::sexp r_dt, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_seed, cpp11::sexp r_deterministic) {
  return dust2::r::dust2_cpu_alloc<{{class}}>(r_pars, r_time, r_dt, r_n_particles, r_n_groups, r_seed, r_deterministic);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_run_steps(cpp11::sexp ptr, cpp11::sexp r_n_steps) {
  return dust2::r::dust2_cpu_run_steps<{{class}}>(ptr, r_n_steps);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_cpu_run_to_time<{{class}}>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_state(cpp11::sexp ptr, bool grouped) {
  return dust2::r::dust2_cpu_state<{{class}}>(ptr, grouped);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_time(cpp11::sexp ptr) {
  return dust2::r::dust2_cpu_time<{{class}}>(ptr);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_set_state_initial(cpp11::sexp ptr) {
  return dust2::r::dust2_cpu_set_state_initial<{{class}}>(ptr);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_set_state(cpp11::sexp ptr, cpp11::sexp r_state, bool grouped) {
  return dust2::r::dust2_cpu_set_state<{{class}}>(ptr, r_state, grouped);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  return dust2::r::dust2_cpu_reorder<{{class}}>(ptr, r_index);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_cpu_rng_state<{{class}}>(ptr);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  return dust2::r::dust2_cpu_set_rng_state<{{class}}>(ptr, r_rng_state);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_cpu_set_time<{{class}}>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_update_pars(cpp11::sexp ptr, cpp11::list pars, bool grouped) {
  return dust2::r::dust2_cpu_update_pars<{{class}}>(ptr, pars, grouped);
}

[[cpp11::register]]
SEXP dust2_cpu_{{name}}_simulate(cpp11::sexp ptr, cpp11::sexp r_times, cpp11::sexp r_index, bool grouped) {
  return dust2::r::dust2_cpu_simulate<{{class}}>(ptr, r_times, r_index, grouped);
}
