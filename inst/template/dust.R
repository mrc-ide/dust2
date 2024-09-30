{{name}} <- function() {
  package <- list(name = "{{package}}", path = "{{{{path_pkg}}}}")
  dust2::dust_system_generator("{{name}}", "{{time_type_property}}", {{default_dt}},
                               package)
}
