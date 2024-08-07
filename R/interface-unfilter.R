##' Create an "unfilter" object, which can be used to compute a
##' deterministic likelihood following the same algorithm as the
##' particle filter, but limited to a single particle.  The name for
##' this method will change in future.
##'
##' @title Create an unfilter
##'
##' @inheritParams dust_system_create
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
dust_unfilter_create <- function(generator, time_start, data,
                                 n_particles = 1, n_groups = NULL,
                                 dt = 1, n_threads = 1, index_state = NULL,
                                 preserve_particle_dimension = FALSE,
                                 preserve_group_dimension = FALSE) {
  call <- environment()
  check_generator_for_filter(generator, "unfilter", call = call)
  assert_scalar_size(n_particles, allow_zero = FALSE, call = call)
  assert_scalar_logical(preserve_particle_dimension, call = call)

  data <- prepare_data(data, n_groups, call = call)
  time_start <- check_time_start(time_start, data$time, call = call)
  dt <- check_dt(dt, call = call)

  n_groups <- data$n_groups
  preserve_group_dimension <- preserve_group_dimension || n_groups > 1
  preserve_particle_dimension <- preserve_particle_dimension || n_particles > 1

  index_state <- check_index(index_state, call = call)
  n_threads <- check_n_threads(n_threads, n_particles, n_groups)

  inputs <- list(time_start = time_start,
                 time = data$time,
                 dt = dt,
                 data = data$data,
                 n_particles = n_particles,
                 n_groups = n_groups,
                 n_threads = n_threads,
                 index_state = index_state,
                 preserve_particle_dimension = preserve_particle_dimension,
                 preserve_group_dimension = preserve_group_dimension)

  res <- list2env(
    list(inputs = inputs,
         n_particles = as.integer(n_particles),
         n_groups = as.integer(n_groups),
         deterministic = TRUE,
         methods = generator$methods$unfilter,
         index_state = index_state,
         preserve_particle_dimension = preserve_particle_dimension,
         preserve_group_dimension = preserve_group_dimension),
    parent = emptyenv())
  class(res) <- "dust_unfilter"
  res
}


unfilter_create <- function(unfilter, pars) {
  inputs <- unfilter$inputs
  list2env(
    unfilter$methods$alloc(pars,
                           inputs$time_start,
                           inputs$time,
                           inputs$dt,
                           inputs$data,
                           inputs$n_particles,
                           inputs$n_groups,
                           inputs$n_threads,
                           inputs$index_state),
    unfilter)
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
##' @param adjoint Logical, indicating if we should enable adjoint
##'   history saving.  This can be enabled even when your model does not
##'   support adjoints!  But you will not be able to compute
##'   gradients.
##'
##' @inheritParams dust_filter_run
##'
##' @return A vector of likelihood values, with as many elements as
##'   there are groups.
##'
##' @export
dust_unfilter_run <- function(unfilter, pars, initial = NULL,
                              save_history = FALSE, adjoint = FALSE,
                              index_group = NULL) {
  check_is_dust_unfilter(unfilter)
  # TODO: parameter checking and updating need to go with our index too.
  if (!is.null(pars)) {
    pars <- check_pars(pars, unfilter$n_groups,
                       unfilter$preserve_group_dimension)
  }
  if (is.null(unfilter$ptr)) {
    if (is.null(pars)) {
      cli::cli_abort("'pars' cannot be NULL, as unfilter is not initialised",
                     arg = "pars")
    }
    unfilter_create(unfilter, pars)
  } else if (!is.null(pars)) {
    unfilter$methods$update_pars(unfilter$ptr, pars, index_group)
  }
  unfilter$methods$run(unfilter$ptr,
                       initial,
                       save_history,
                       adjoint,
                       index_group,
                       unfilter$preserve_particle_dimension,
                       unfilter$preserve_group_dimension)
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
dust_unfilter_last_history <- function(unfilter, index_group = NULL) {
  check_is_dust_unfilter(unfilter)
  if (is.null(unfilter$ptr)) {
    cli::cli_abort(c(
      "History is not current",
      i = "Unfilter has not yet been run"))
  }
  unfilter$methods$last_history(unfilter$ptr,
                                index_group,
                                unfilter$preserve_particle_dimension,
                                unfilter$preserve_group_dimension)
}


##' Fetch the last gradient created by running an unfilter.  This
##' errors if the last call to [dust_unfilter_run] did not use
##' `adjoint = TRUE`.  The first time you call this (after a
##' particular set of parameters) it will trigger running the reverse
##' model.
##'
##' @title Fetch last unfilter gradient
##'
##' @inheritParams dust_unfilter_run
##'
##' @return A vector (if ungrouped) or a matrix (if grouped).
##'
##' @export
dust_unfilter_last_gradient <- function(unfilter, index_group = NULL) {
  check_is_dust_unfilter(unfilter)
  if (is.null(unfilter$ptr)) {
    cli::cli_abort(c(
      "Gradient is not current",
      i = "Unfilter has not yet been run"))
  }
  unfilter$methods$last_gradient(unfilter$ptr,
                                 index_group,
                                 unfilter$preserve_particle_dimension,
                                 unfilter$preserve_group_dimension)
}


check_is_dust_unfilter <- function(unfilter, call = parent.frame()) {
  if (!inherits(unfilter, "dust_unfilter")) {
    cli::cli_abort("Expected 'unfilter' to be a 'dust_unfilter' object",
                   arg = "unfilter", call = call)
  }
}
