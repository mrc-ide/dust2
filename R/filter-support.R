check_generator_for_filter <- function(generator, what, call = NULL) {
  check_is_dust_system_generator(generator)
  if (!generator$properties$has_compare) {
    cli::cli_abort(
      paste("Can't create {what}; the '{generator$name}' system does",
            "not have 'compare_data' support"),
      arg = "generator")
  }
  generator
}


check_dt <- function(dt, call = NULL) {
  assert_scalar_numeric(dt, call = call)
  if (dt <= 0) {
    cli::cli_abort("Expected 'dt' to be greater than 0")
  }
  if (dt > 1) {
    cli::cli_abort("Expected 'dt' to be at most 1")
  }
  if (!rlang::is_integerish(1 / dt)) {
    cli::cli_abort("Expected 'dt' to be the inverse of an integer",
                   arg = "dt", call = call)
  }
  dt
}


check_index <- function(index, max = NULL, unique = FALSE,
                        name = deparse(substitute(index)), call = NULL) {
  if (!is.null(index)) {
    assert_integer(index, name = name, call = call)
    if (any(index < 1)) {
      cli::cli_abort("All elements of '{name}' must be at least 1",
                     arg = name, call = call)
    }
    if (!is.null(max) && any(index > max)) {
      cli::cli_abort("All elements of '{name}' must be at most {max}",
                     arg = name, call = call)
    }
    if (unique && anyDuplicated(index)) {
      cli::cli_abort("All elements of '{name}' must be distinct",
                     arg = name, call = call)
    }
  }
  index
}


check_time_start <- function(time_start, time, call = NULL) {
  assert_scalar_integer(time_start, call = call)
  if (time_start > time[[1]]) {
    cli::cli_abort(
      paste("'time_start' ({time_start}) is later than the first time",
            "in 'data' ({time[[1]]})"),
      call = call)
  }
  time_start
}
