[[cpp11::register]]
SEXP dust2_system_{{name}}_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_run_to_time<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_state(cpp11::sexp ptr, cpp11::sexp r_index_state, cpp11::sexp r_index_particle, cpp11::sexp r_index_group, bool preserve_group_dimension) {
  return dust2::r::dust2_system_state<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_index_state, r_index_particle, r_index_group, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_time(cpp11::sexp ptr) {
  return dust2::r::dust2_system_time<dust2::dust_{{time_type}}<{{class}}>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_set_state_initial(cpp11::sexp ptr) {
  return dust2::r::dust2_system_set_state_initial<dust2::dust_{{time_type}}<{{class}}>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_set_state(cpp11::sexp ptr, cpp11::sexp r_state, bool preserve_group_dimension) {
  return dust2::r::dust2_system_set_state<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_state, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  return dust2::r::dust2_system_reorder<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_index);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_system_rng_state<dust2::dust_{{time_type}}<{{class}}>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  return dust2::r::dust2_system_set_rng_state<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_rng_state);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_set_time<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_update_pars(cpp11::sexp ptr, cpp11::list pars) {
  return dust2::r::dust2_system_update_pars<dust2::dust_{{time_type}}<{{class}}>>(ptr, pars);
}

[[cpp11::register]]
SEXP dust2_system_{{name}}_simulate(cpp11::sexp ptr, cpp11::sexp r_times, cpp11::sexp r_index, bool preserve_group_dimension) {
  return dust2::r::dust2_system_simulate<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_times, r_index, preserve_group_dimension);
}
