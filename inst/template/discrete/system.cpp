[[cpp11::register]]
SEXP dust2_system_{{name}}_alloc(cpp11::list r_pars, cpp11::sexp r_time, cpp11::sexp r_dt, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_seed, cpp11::sexp r_deterministic, cpp11::sexp r_n_threads) {
  return dust2::r::dust2_discrete_alloc<{{class}}>(r_pars, r_time, r_dt, r_n_particles, r_n_groups, r_seed, r_deterministic, r_n_threads);
}
