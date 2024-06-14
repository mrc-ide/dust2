[[cpp11::register]]
SEXP dust2_system_{{name}}_compare_data(cpp11::sexp ptr, cpp11::sexp r_data, bool grouped) {
  return dust2::r::dust2_system_compare_data<dust2::dust_{{time_type}}<{{class}}>>(ptr, r_data, grouped);
}
