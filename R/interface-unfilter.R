##' Create an "unfilter" object, which can be used to compute a
##' deterministic likelihood following the same algorithm as the
##' particle filter, but limited to a single particle.  The name for
##' this method will change in future.
##'
##' @title Create an unfilter
##'
##' @inheritParams dust_filter_create
##'
##' @param n_particles The number of particles to run.  Typically this
##'   is 1, but you can run with more than 1 if you want - currently
##'   they produce the same likelihood but if you provide different
##'   initial conditions then you would see different likelihoods.
##'
##' @return A `dust_unfilter` object, which can be used with
##'   [dust_unfilter_run]
##'
##' @export
dust_unfilter_create <- function(generator, time_start, time, data,
                                 n_particles = 1, n_groups = 0,
                                 dt = 1, index = NULL) {
  call <- environment()
  check_generator_for_filter(generator, "unfilter", call = call)
  assert_scalar_size(n_particles, call = call)
  assert_scalar_size(n_groups, allow_zero = TRUE, call = call)
  check_time_sequence(time_start, time, call = call)
  check_dt(dt, call = call)
  check_data(data, length(time), n_groups, call = call)
  check_index(index, call = call)

  create <- function(unfilter, pars) {
    list2env(
      generator$methods$unfilter$alloc(pars, time_start, time, dt, data,
                                       n_particles, n_groups, index),
      unfilter)
    unfilter$create <- NULL
  }
  res <- list2env(
    list(create = create,
         n_particles = as.integer(n_particles),
         n_groups = as.integer(max(n_groups), 1),
         deterministic = TRUE,
         methods = generator$methods$unfilter,
         index = index),
    parent = emptyenv())
  class(res) <- "dust_unfilter"
  res
}


##' Run unfilter
##'
##' @title Run unfilter
##'
##' @inheritParams dust_filter_run
##'
##' @param unfilter A `dust_unfilter` object, created by
##'   [dust_unfilter_create]
##'
##' @inheritParams dust_filter_run
##'
##' @return A vector of likelihood values, with as many elements as
##'   there are groups.
##'
##' @export
dust_unfilter_run <- function(unfilter, pars, initial = NULL,
                              save_history = FALSE) {
  check_is_dust_unfilter(unfilter)
  if (is.null(unfilter$ptr)) {
    if (is.null(pars)) {
      cli::cli_abort("'pars' cannot be NULL, as unfilter is not initialised",
                     arg = "pars")
    }
    unfilter$create(unfilter, pars)
  } else if (!is.null(pars)) {
    unfilter$methods$update_pars(unfilter$ptr, pars, unfilter$grouped)
  }
  unfilter$methods$run(unfilter$ptr, initial, save_history, unfilter$grouped)
}


##' Fetch the last history created by running an unfilter.  This
##' errors if the last call to [dust_unfilter_run] did not use
##' `save_history = TRUE`.
##'
##' @title Fetch last unfilter history
##'
##' @inheritParams dust_unfilter_run
##'
##' @return An array
##'
##' @export
dust_unfilter_last_history <- function(unfilter) {
  check_is_dust_unfilter(unfilter)
  if (is.null(unfilter$ptr)) {
    cli::cli_abort(c(
      "History is not current",
      i = "Unfilter has not yet been run"))
  }
  unfilter$methods$last_history(unfilter$ptr, unfilter$grouped)
}


check_is_dust_unfilter <- function(unfilter, call = parent.frame()) {
  if (!inherits(unfilter, "dust_unfilter")) {
    cli::cli_abort("Expected 'unfilter' to be a 'dust_unfilter' object",
                   arg = "unfilter", call = call)
  }
}
