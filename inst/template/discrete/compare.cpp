[[cpp11::register]]
SEXP dust2_unfilter_{{name}}_alloc(cpp11::list r_pars, cpp11::sexp r_time_start, cpp11::sexp r_time, cpp11::sexp r_dt, cpp11::list r_data, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_index) {
  return dust2::r::dust2_discrete_unfilter_alloc<{{class}}>(r_pars, r_time_start, r_time, r_dt, r_data, r_n_particles, r_n_groups, r_index);
}

[[cpp11::register]]
SEXP dust2_filter_{{name}}_alloc(cpp11::list r_pars, cpp11::sexp r_time_start, cpp11::sexp r_time, cpp11::sexp r_dt, cpp11::list r_data, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_index, cpp11::sexp r_seed) {
  return dust2::r::dust2_discrete_filter_alloc<{{class}}>(r_pars, r_time_start, r_time, r_dt, r_data, r_n_particles, r_n_groups, r_index, r_seed);
}
