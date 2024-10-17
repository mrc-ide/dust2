[[cpp11::register]]
SEXP dust2_unfilter_{{name}}_last_gradient(cpp11::sexp ptr, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_discrete_unfilter_last_gradient<dust2::dust_{{time_type}}<{{class}}>>(ptr, preserve_particle_dimension, preserve_group_dimension);
}
