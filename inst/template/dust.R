{{name}} <- function() {
  if (isNamespaceLoaded("{{package}}")) {
    env <- asNamespace("{{package}}")
  } else {
    env <- dust2:::load_temporary_package("{{{{path_pkg}}}}", "{{package}}", FALSE)
  }
  dust2::dust_system_generator("{{name}}", "{{time_type_property}}", {{default_dt}}, env)
}
