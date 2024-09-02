##' Prepare data for use with [dust_unfilter_create] or
##' [dust_filter_create].  You do not have to use this function if you
##' name your data.frame with our standard column names (i.e., `time`
##' column containing the time) as it will be called within the filter
##' functions directly.  However, you can use this to validate your
##' data separately or to use different columns than the defaults.
##'
##' @title Prepare data
##'
##' @param data A data.frame containing time and data to fit.  By
##'   default we expect a column `time` (or one with the name given as
##'   the argument `time`) and one or more columns of data to fit to.
##'
##' @param time Optional name of a column within `data` to use for
##'   time.
##'
##' @param group Optional name of a column within `data` to use for
##'   groups
##'
##' @return A data.frame, with the addition of the class attribute
##'   `dust_filter_data`; once created you should not modify this
##'   object.
##'
##' @export
dust_filter_data <- function(data, time = NULL, group = NULL) {
  if (inherits(data, "dust_filter_data")) {
    return(data)
  }

  assert_is(data, "data.frame")
  if (nrow(data) == 0) {
    cli::cli_abort("Expected 'data' to have at least one row")
  }

  ## First check time; this has to be nailed down right away:
  time <- time %||% "time"
  if (is.null(data[[time]])) {
    cli::cli_abort(
      "Did not find column '{time}' in 'data'")
  }
  ## TODO: we should support noninteger time here but that will take a
  ## bit of a fiddle still...
  assert_integer(data[[time]], name = sprintf('data[["%s"]]', time),
                 arg = "data")

  is_plausibly_grouped <- anyDuplicated(data[[time]]) > 0
  if (is_plausibly_grouped && is.null(group)) {
    if (is.null(data[["group"]])) {
      cli::cli_abort(
        c(paste("'data' looks grouped, but 'group' not specified and column",
                "'group' not found"),
          i = paste("You have duplicated times in 'data[[\"{time}\"]]',",
                    "so I believe 'data' is grouped but you have not provided",
                    "the 'group' argument")))
    }
    group <- "group"
  } else if (!is.null(group)) {
    assert_scalar_character(group)
    if (is.null(data[[group]])) {
      cli::cli_abort(
        c("Did not find column '{group}' in 'data'",
          i = "You provided the argument 'group'"))
    }
  }

  if (is.null(group)) {
    n_groups <- 1
    data <- data[order(data[[time]]), ]
  } else {
    assert_integer(data[[group]], name = sprintf('data[["%s"]]', group),
                   arg = "data")
    groups <- unique(data[[group]])
    n_groups <- length(groups)
    if (!setequal(groups, seq_len(n_groups))) {
      cli::cli_abort(
        paste("Expected 'data[[\"{group}\"]]' to contain integer values in",
              "[1, {n_groups}], and contain all those values"))
    }

    ## Detect balance here; this is not super easy, and in particular
    ## not easy to report back errors about. Let's refine the error
    ## reporting later.
    err <- length(unique(lapply(split(data[[time]], data[[group]]), sort))) > 1
    if (err) {
      cli::cli_abort(
        c("Not all groups in 'data' have the same times",
          i = paste("Each group (from the '{group}' column) must have",
                    "the same times (from the '{time}' column), but",
                    "your groups have different times")))
    }
    data <- data[order(data[[time]], data[[group]]), ]
  }

  if (length(setdiff(names(data), c(time, group))) == 0) {
    cli::cli_abort(
      paste("Expected 'data' to have at least one column in addition to",
            "{squote(c(time, group))}"))
  }

  rownames(data) <- NULL
  attr(data, "time") <- time
  attr(data, "group") <- group
  attr(data, "n_groups") <- n_groups
  class(data) <- c("dust_filter_data", class(data))
  data
}


prepare_data <- function(data, n_groups, call = NULL) {
  if (!is.null(n_groups)) {
    assert_scalar_integer(n_groups, call = call)
  }
  data <- dust_filter_data(data)

  name_time <- attr(data, "time")
  name_group <- attr(data, "group")
  n_groups_data <- attr(data, "n_groups")

  if (is.null(n_groups)) {
    n_groups <- n_groups_data
  } else if (n_groups != n_groups_data) {
    cli::cli_abort(
      "Expected 'data' to have {n_groups} groups, but it had {n_groups_data}",
      call = call)
  }

  time <- data[[name_time]]
  if (!is.null(name_group)) {
    time <- time[data[[name_group]] == 1]
  }

  cols <- names(data) != name_time
  data_by_time <- unname(split(data[cols], data[[name_time]]))
  is_list_column <- vlapply(data[cols], is.list)
  row_as_list <- function(x) {
    ret <- as.list(x)
    if (any(is_list_column)) {
      stopifnot(all(lengths(x[is_list_column]) == 1))
      ret[is_list_column] <- lapply(ret[is_list_column], "[[", 1)
    }
    ret
  }

  if (is.null(name_group)) {
    data_by_time <- lapply(data_by_time, function(x) list(row_as_list(x)))
  } else {
    data_by_time <- lapply(data_by_time, function(el) {
      lapply(unname(
        split(el[names(el) != name_group], el[[name_group]])),
        row_as_list)
    })
  }

  list(time = time,
       n_groups = n_groups_data,
       data = data_by_time)
}
