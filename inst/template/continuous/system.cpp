[[cpp11::register]]
SEXP dust2_system_{{name}}_alloc(cpp11::list r_pars, cpp11::sexp r_time, cpp11::list r_time_control, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_seed, cpp11::sexp r_deterministic, cpp11::sexp r_n_threads) {
  return dust2::r::dust2_continuous_alloc<{{class}}>(r_pars, r_time, r_time_control, r_n_particles, r_n_groups, r_seed, r_deterministic, r_n_threads);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_internals(cpp11::sexp ptr, bool include_coefficients) {
  return dust2::r::dust2_system_internals<dust2::dust_{{time_type}}<{{class}}>>(ptr, include_coefficients);
}
