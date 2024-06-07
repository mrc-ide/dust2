[[cpp11::register]]
SEXP dust2_discrete_{{name}}_unfilter_last_gradient(cpp11::sexp ptr, bool grouped) {
  return dust2::r::dust2_discrete_unfilter_last_gradient<{{class}}>(ptr, grouped);
}
