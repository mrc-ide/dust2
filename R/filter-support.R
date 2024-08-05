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


check_time_sequence <- function(time_start, time, call = NULL) {
  assert_scalar_integer(time_start, call = call)
  assert_integer(time, call = call)
  if (length(time) == 0) {
    cli::cli_abort("Expected at least one value in 'time'",
                   arg = "time", call = call)
  }
  err <- diff(c(time_start, time)) <= 0
  if (any(err)) {
    i <- which(err)
    nm_other <- ifelse(i == 1, "time_start", sprintf("time[%d]", i - 1))
    detail <- sprintf("'time[%d]' (%d) must be greater than '%s' (%d)",
                      i, time[i], nm_other, c(time_start, time)[i])
    if (length(detail) > 5) {
      detail <- c(
        detail[1:4],
        sprintf("...and %d other errors", sum(err) - 4))
    }
    cli::cli_abort(
      c("Time sequence is not strictly increasing",
        set_names(detail, "x")),
      arg = "time", call = call)
  }
  as.numeric(time)
}


check_data <- function(data, n_time, n_groups, preserve_group_dimension,
                       call = NULL) {
  assert_list(data, call = call)
  assert_length(data, n_time, call = call)
  if (preserve_group_dimension) {
    len <- lengths(data)
    err <- len != n_groups
    if (any(err)) {
      detail <- sprintf("Error for element %d, which has length %d",
                        which(err), len[err])
      if (length(detail) > 5) {
        detail <- c(
          detail[1:4],
          sprintf("...and %d other elements", sum(err) - 4))
      }
      if (n_groups > 1) {
        justification <- "'n_groups' is greater than one"
      } else {
        justification <- "preserve_group_dimension was TRUE"
      }
      cli::cli_abort(
        c("Expected all elements of 'data' to have length {n_groups}",
          i = paste(
            "You have a grouped system ({justification})",
            "so each element in data must be a list with data for each group",
            "in turn"),
          set_names(detail, "x")),
        arg = "data", call = call)
    }
  } else {
    data <- lapply(data, function(el) list(el))
  }
  data
}


check_index <- function(index, call = NULL) {
  if (!is.null(index)) {
    assert_integer(index, call = call)
    if (any(index < 0)) {
      cli::cli_abort("All elements of 'index' must be at least 1",
                     arg = "index", call = call)
    }
  }
  index
}
