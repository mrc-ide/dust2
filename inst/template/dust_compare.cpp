[[cpp11::register]]
SEXP dust2_discrete_{{name}}_compare_data(cpp11::sexp ptr, cpp11::sexp r_data, bool grouped) {
  return dust2::r::dust2_discrete_compare_data<{{class}}>(ptr, r_data, grouped);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_unfilter_alloc(cpp11::list r_pars, cpp11::sexp r_time_start, cpp11::sexp r_time, cpp11::sexp r_dt, cpp11::list r_data, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_index) {
  return dust2::r::dust2_discrete_unfilter_alloc<{{class}}>(r_pars, r_time_start, r_time, r_dt, r_data, r_n_particles, r_n_groups, r_index);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_unfilter_update_pars(cpp11::sexp ptr, cpp11::list r_pars, bool grouped) {
  return dust2::r::dust2_discrete_unfilter_update_pars<{{class}}>(ptr, r_pars, grouped);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_unfilter_run(cpp11::sexp ptr, cpp11::sexp r_initial, bool save_history, bool adjoint, bool grouped) {
  return dust2::r::dust2_discrete_unfilter_run<{{class}}>(ptr, r_initial, save_history, adjoint, grouped);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_unfilter_last_history(cpp11::sexp ptr, bool grouped) {
  return dust2::r::dust2_discrete_unfilter_last_history<{{class}}>(ptr, grouped);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_filter_alloc(cpp11::list r_pars, cpp11::sexp r_time_start, cpp11::sexp r_time, cpp11::sexp r_dt, cpp11::list r_data, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_index, cpp11::sexp r_seed) {
  return dust2::r::dust2_discrete_filter_alloc<{{class}}>(r_pars, r_time_start, r_time, r_dt, r_data, r_n_particles, r_n_groups, r_index, r_seed);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_filter_update_pars(cpp11::sexp ptr, cpp11::list r_pars, bool grouped) {
  return dust2::r::dust2_discrete_filter_update_pars<{{class}}>(ptr, r_pars, grouped);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_filter_run(cpp11::sexp ptr, cpp11::sexp r_initial, bool save_history, bool grouped) {
  return dust2::r::dust2_discrete_filter_run<{{class}}>(ptr, r_initial, save_history, grouped);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_filter_last_history(cpp11::sexp ptr, bool grouped) {
  return dust2::r::dust2_discrete_filter_last_history<{{class}}>(ptr, grouped);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_filter_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_discrete_filter_rng_state<{{class}}>(ptr);
}

[[cpp11::register]]
SEXP dust2_discrete_{{name}}_filter_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  return dust2::r::dust2_discrete_filter_set_rng_state<{{class}}>(ptr, r_rng_state);
}
