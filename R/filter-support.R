check_generator_for_filter <- function(generator, what, call = NULL) {
  check_is_dust_system_generator(generator, substitute(generator))
  if (!attr(generator, "properties")$has_compare) {
    name <- attr(generator, "name")
    cli::cli_abort(
      paste("Can't create {what}; the '{name}' system does",
            "not have 'compare_data' support"),
      arg = "generator")
  }
  generator
}


check_system_dt <- function(dt, generator, name = "dt", call = NULL) {
  time_type <- properties <- attr(generator, "properties")$time_type
  if (time_type == "continuous") {
    if (!is.null(dt)) {
      cli::cli_abort("Can't use '{name}' with continuous-time systems",
                     call = call)
    }
  } else {
    check_dt(dt %||% attr(generator, "default_dt"),
             allow_null = time_type == "mixed",
             name = name,
             call = call)
  }
}


check_system_ode_control <- function(ode_control, generator,
                                     name = "ode_control", call = NULL) {
  time_type <- attr(generator, "properties")$time_type
  if (time_type == "discrete") {
    if (!is.null(ode_control)) {
      cli::cli_abort("Can't use 'ode_control' with discrete-time systems")
    }
  } else {
    if (is.null(ode_control)) {
      ode_control <- dust_ode_control()
    } else {
      assert_is(ode_control, "dust_ode_control", call = environment())
    }
  }
}


check_dt <- function(dt, allow_null = FALSE, name = deparse(substitute(dt)),
                     call = NULL) {
  if (allow_null && is.null(dt)) {
    return(dt)
  }
  assert_scalar_numeric(dt, call = call)
  if (dt <= 0) {
    cli::cli_abort("Expected '{name}' to be greater than 0")
  }
  if (dt > 1) {
    cli::cli_abort("Expected '{name}' to be at most 1")
  }
  if (!rlang::is_integerish(1 / dt)) {
    cli::cli_abort("Expected '{name}' to be the inverse of an integer",
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
  ## TODO: port over bits from check_time() here...
  assert_scalar_integer(time_start, call = call)
  if (time_start > time[[1]]) {
    cli::cli_abort(
      paste("'time_start' ({time_start}) is later than the first time",
            "in 'data' ({time[[1]]})"),
      call = call)
  }
  time_start
}
