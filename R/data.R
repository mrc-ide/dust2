##' Collect time and data inputs for running a particle filter
##' ([dust_filter_create]) or unfilter ([dust_unfilter_create]).  This
##' function collects related elements that need validation.
##'
##' @title Create dust data
##'
##' @param time_start The start time for the simulation - this is
##'   typically before the first data point.  Must be an integer-like
##'   value.
##'
##' @param time A vector of times, each of which has a corresponding
##'   entry in `data`.  The model will stop at each of these times to
##'   compute the likelihood using the compare function.
##'
##' @param data The data to compare against.  This must be a list with
##'   the same length as `time`, each element of which corresponds to
##'   the data required for the model.  If the model is ungrouped then
##'   each element of `data` is a list with elements corresponding to
##'   whatever your model requires.  If your model is grouped, this
##'   should be a list with as many elements as your model has groups,
##'   with each element corresponding to the data your model requires.
##'   We will likely introduce a friendlier data.frame based input
##'   soon.
##'
##' @return An object of class `dust_data`, which should not be
##'   modified after creation.  This can be passed as `data` to
##'   `dust_filter_create` and `dust_unfilter_create`.
##' 
##' @export
dust_data <- function(time_start, time, data, dt = 1, n_groups = 0) {
  call <- rlang::current_call()

  assert_scalar_integer(time_start, call = call) # check_time
  assert_scalar_numeric(dt, call = call)
  if (!rlang::is_integerish(1 / dt)) {
    cli::cli_abort("Expected 'dt' to be the inverse of an integer",
                   arg = "dt", call = call)
  }

  assert_integer(time)
  if (length(time) == 0) {
    cli::cli_abort("Expected 'time' to have at least one value",
                   arg = "time", call = call)
  }

  time_all <- c(time_start, time)
  err <- diff(time_all) <= 0
  if (any(err)) {
    i <- which(err)[[1]]
    curr <- time_all[[i + 1]]
    prev <- time_all[[i]]
    hint <- NULL
    if (i == 1) {
      hint <- c(i = "Previous value comes from 'time_start'")
    }
    if (sum(err) > 1) {
      hint <- c(x = "And {sum(i) - 1} other invalid time{?s}")
    }
    cli::cli_abort(
      c(paste("Expected 'time[{i}]' ({curr}) to be larger than",
              "the previous value ({prev})"),
        hint),
      call = "time", call = call)
  }

  ## Might also support a list-matrix here too?
  if (is.data.frame(data)) {
    cli::cli_abort("A data.frame for 'data' is not yet supported",
                   arg = "data", call = call)
  }
  if (length(data) != length(time)) {
    cli::cli_abort(
      c("Expected 'data' to have length {length(time)}",
        i = "'data' must be a list with the same length as 'time'"),
      arg = "data", call = call)
  }
  grouped <- n_groups > 0
  if (grouped) {
    len <- lengths(data)
    err <- len != n_groups
    if (any(err)) {
      if (all(err)) {
        detail <- "Error for all {length(data)} entries"
      } else if (sum(err) > 0) {
        detail <- "Error for elements {which(i)}"
      }
      cli::cli_abort(
        c("All elements of 'data' must have '{n_groups}' elements",
          x = detail,
          i = "The number of elements at each time must match 'n_groups'"),
        arg = "data", call = call)
    }
  }

  ret <- list(time_start = time_start,
              time = time,
              dt = dt,
              data = data,
              n_groups = n_groups,
              grouped = grouped)
  clas(ret) <- "dust_data"
  ret
}
