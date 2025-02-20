[[cpp11::register]]
SEXP dust2_system_{{name}}_compare_data(cpp11::sexp ptr, cpp11::list r_data, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_system_compare_data<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_data, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_unfilter_{{name}}_update_pars(cpp11::sexp ptr, cpp11::list r_pars, cpp11::sexp r_index_group) {
  return dust2::r::dust2_unfilter_update_pars<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_pars, r_index_group);
}

[[cpp11::register]]
SEXP dust2_unfilter_{{name}}_run(cpp11::sexp ptr, cpp11::sexp r_initial, bool save_trajectories, cpp11::sexp save_snapshots, bool adjoint, cpp11::sexp r_index_state, cpp11::sexp r_index_group, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_unfilter_run<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_initial, save_trajectories, save_snapshots, adjoint, r_index_state, r_index_group, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_unfilter_{{name}}_last_trajectories(cpp11::sexp ptr, bool select_random_particle, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_unfilter_last_trajectories<dust2::dust_{{time_type}}<{{class}}>>(ptr, select_random_particle, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_unfilter_{{name}}_last_snapshots(cpp11::sexp ptr, bool select_random_particle, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_unfilter_last_snapshots<dust2::dust_{{time_type}}<{{class}}>>(ptr, select_random_particle, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_unfilter_{{name}}_last_state(cpp11::sexp ptr, bool select_random_particle, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_unfilter_last_state<dust2::dust_{{time_type}}<{{class}}>>(ptr, select_random_particle, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_filter_{{name}}_update_pars(cpp11::sexp ptr, cpp11::list r_pars, cpp11::sexp r_index_group) {
  return dust2::r::dust2_filter_update_pars<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_pars, r_index_group);
}

[[cpp11::register]]
SEXP dust2_filter_{{name}}_run(cpp11::sexp ptr, cpp11::sexp r_initial, bool save_trajectories, cpp11::sexp save_snapshots, bool adjoint, cpp11::sexp index_state, cpp11::sexp index_group, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_filter_run<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_initial, save_trajectories, save_snapshots, adjoint, index_state, index_group, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_filter_{{name}}_last_trajectories(cpp11::sexp ptr, bool select_random_particle, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_filter_last_trajectories<dust2::dust_{{time_type}}<{{class}}>>(ptr, select_random_particle, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_filter_{{name}}_last_snapshots(cpp11::sexp ptr, bool select_random_particle, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_filter_last_snapshots<dust2::dust_{{time_type}}<{{class}}>>(ptr, select_random_particle, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_filter_{{name}}_last_state(cpp11::sexp ptr, bool select_random_particle, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_filter_last_state<dust2::dust_{{time_type}}<{{class}}>>(ptr, select_random_particle, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_filter_{{name}}_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_filter_rng_state<dust2::dust_{{time_type}}<{{class}}>>(ptr);
}

[[cpp11::register]]
SEXP dust2_filter_{{name}}_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  return dust2::r::dust2_filter_set_rng_state<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_rng_state);
}
