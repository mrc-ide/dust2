[[cpp11::register]]
SEXP dust2_system_{{name}}_alloc(cpp11::list r_pars, cpp11::sexp r_time, cpp11::sexp r_control, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_seed, cpp11::sexp r_deterministic) {
  return dust2::r::dust2_continuous_alloc<{{class}}>(r_pars, r_time, r_control, r_n_particles, r_n_groups, r_seed, r_deterministic);
}
