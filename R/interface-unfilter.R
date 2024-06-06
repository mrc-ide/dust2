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
dust_unfilter_create <- function(generator, pars, time_start, time, data,
                                 n_particles = 1, n_groups = 0,
                                 dt = 1, index = NULL) {
  check_is_dust_system_generator(generator)
  if (!generator$properties$has_compare) {
    ## This moves into something general soon?
    cli::cli_abort(
      paste("Can't create unfilter; the '{generator$name}' system does",
            "not have 'compare_data' support"),
      arg = "generator")
  }
  res <- generator$methods$unfilter$alloc(pars, time_start, time, dt, data,
                                          n_particles, n_groups, index)
  res$name <- generator$name
  res$n_particles <- as.integer(n_particles)
  res$n_groups <- as.integer(max(n_groups), 1)
  res$deterministic <- TRUE
  res$methods <- generator$methods$unfilter
  res$index <- index
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
dust_unfilter_run <- function(unfilter, pars = NULL, initial = NULL,
                              save_history = FALSE) {
  check_is_dust_unfilter(unfilter)
  unfilter$methods$run(unfilter$ptr, pars, initial, save_history,
                       unfilter$grouped)
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
  unfilter$methods$last_history(unfilter$ptr, unfilter$grouped)
}


check_is_dust_unfilter <- function(unfilter, call = parent.frame()) {
  if (!inherits(unfilter, "dust_unfilter")) {
    cli::cli_abort("Expected 'unfilter' to be a 'dust_unfilter' object",
                   arg = "unfilter", call = call)
  }
}
