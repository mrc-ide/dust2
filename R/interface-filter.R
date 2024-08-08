##' Create a particle filter object
##'
##' @title Create a particle filter
##'
##' @param generator A system generator object, with class
##'   `dust_system_generator`.  The system must support `compare_data`
##'   to be used with this function.
##'
##' @param time_start The start time for the simulation - this is
##'   typically before the first data point.  Must be an integer-like
##'   value.
##'
##' @param data The data to fit to.  This can be a `data.frame`, in
##'   which case it will be passed into [dust_filter_data] for
##'   validation, or it can be a [dust_filter_data]-augmented
##'   data.frame.  The times for comparison will be taken from this,
##'   and `time_start` must be earlier than the earliest time.
##'
##' @param n_particles The number of particles to run.  Larger numbers
##'   will give lower variance in the likelihood estimate but run more
##'   slowly.
##'
##' @param n_groups The number of parameter groups.  If `NULL`, this
##'   will be taken from `data`.  If given, then the number of groups
##'   in `data` will be checked against this number.
##'
##' @inheritParams dust_system_create
##' @inheritParams dust_system_simulate
##'
##' @return A `dust_unfilter` object, which can be used with
##'   [dust_unfilter_run]
##'
##' @export
dust_filter_create <- function(generator, time_start, data,
                               n_particles, n_groups = NULL, dt = 1,
                               index_state = NULL, n_threads = 1,
                               preserve_group_dimension = FALSE,
                               seed = NULL) {
  call <- environment()
  check_generator_for_filter(generator, "filter", call = call)
  assert_scalar_size(n_particles, allow_zero = FALSE, call = call)
  assert_scalar_logical(preserve_group_dimension, call = call)

  data <- prepare_data(data, n_groups, call = call)
  time_start <- check_time_start(time_start, data$time, call = call)
  dt <- check_dt(dt, call = call)

  ## NOTE: there is no preserve_particle_dimension option here because
  ## we will always preserve this dimension.
  n_groups <- data$n_groups
  preserve_group_dimension <- preserve_group_dimension || n_groups > 1

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
                 preserve_group_dimension = preserve_group_dimension)

  res <- list2env(
    list(inputs = inputs,
         initial_rng_state = filter_rng_state(n_particles, n_groups, seed),
         n_particles = n_particles,
         n_groups = n_groups,
         deterministic = FALSE,
         methods = generator$methods$filter,
         index_state = index_state,
         preserve_group_dimension = preserve_group_dimension),
    parent = emptyenv())
  class(res) <- "dust_filter"
  res
}


##' Create an independent copy of a filter.  The new filter is
##' decoupled from the random number streams of the parent filter.  It
##' is also decoupled from the *state size* of the parent filter, so
##' you can use this to create a new filter where the system is
##' fundamentally different but everything else is the same.
##'
##' @title Create copy of filter
##'
##' @inheritParams dust_filter_run
##'
##' @param seed The seed for the filter (see [dust_filter_create])
##'
##' @return A new `dust_filter` object
dust_filter_copy <- function(filter, seed = NULL) {
  dst <- new.env(parent = emptyenv())
  nms <- c("inputs", "n_particles", "n_groups", "deterministic", "methods",
           "index_state", "preserve_group_dimension")
  for (nm in nms) {
    dst[[nm]] <- filter[[nm]]
  }
  dst$initial_rng_state <-
    filter_rng_state(filter$n_particles, filter$n_groups, seed)
  class(dst) <- "dust_filter"
  dst
}


filter_create <- function(filter, pars) {
  inputs <- filter$inputs
  list2env(
    filter$methods$alloc(pars,
                         inputs$time_start,
                         inputs$time,
                         inputs$dt,
                         inputs$data,
                         inputs$n_particles,
                         inputs$n_groups,
                         inputs$n_threads,
                         inputs$index_state,
                         filter$initial_rng_state),
    filter)
  filter$initial_rng_state <- NULL
}


##' Run particle filter
##'
##' @title Run particle filter
##'
##' @param filter A `dust_filter` object, created by
##'   [dust_filter_create]
##'
##' @param pars Optional parameters to run the filter with.  If not
##'   provided, parameters are not updated
##'
##' @param initial Optional initial conditions, as a matrix (state x
##'   particle) or 3d array (state x particle x group).  If not
##'   provided, the system initial conditions are used.
##'
##' @param save_history Logical, indicating if the simulation history
##'   should be saved while the simulation runs; this has a small
##'   overhead in runtime and in memory.  History (particle
##'   trajectories) will be saved at each time in the filter.  If the
##'   filter was constructed using a non-`NULL` `index_state` parameter,
##'   the history is restricted to these states.
##'
##' @return A vector of likelihood values, with as many elements as
##'   there are groups.
##'
##' @export
dust_filter_run <- function(filter, pars, initial = NULL,
                            save_history = FALSE) {
  check_is_dust_filter(filter)
  if (!is.null(pars)) {
    pars <- check_pars(pars, filter$n_groups,
                       filter$preserve_group_dimension)
  }
  if (is.null(filter$ptr)) {
    if (is.null(pars)) {
      cli::cli_abort("'pars' cannot be NULL, as filter is not initialised",
                     arg = "pars")
    }
    filter_create(filter, pars)
  } else if (!is.null(pars)) {
    filter$methods$update_pars(filter$ptr, pars)
  }
  filter$methods$run(filter$ptr, initial, save_history,
                     filter$preserve_group_dimension)
}


##' Fetch the last history created by running a filter.  This
##' errors if the last call to [dust_filter_run] did not use
##' `save_history = TRUE`.
##'
##' @title Fetch last filter history
##'
##' @inheritParams dust_filter_run
##'
##' @return An array
##'
##' @export
dust_filter_last_history <- function(filter) {
  check_is_dust_filter(filter)
  if (is.null(filter$ptr)) {
    cli::cli_abort(c(
      "History is not current",
      i = "Filter has not yet been run"))
  }
  filter$methods$last_history(filter$ptr, filter$preserve_group_dimension)
}


##' Get random number generator (RNG) state from the particle filter.
##'
##' @title Get filter RNG state
##'
##' @inheritParams dust_filter_run
##'
##' @return A raw vector, this could be quite long.  Later we will
##'   describe how you might reseed a filter or system with this state.
##'
##' @export
dust_filter_rng_state <- function(filter) {
  check_is_dust_filter(filter)
  if (is.null(filter$ptr)) {
    filter$initial_rng_state
  } else {
    filter$methods$rng_state(filter$ptr)
  }
}


##' @param rng_state A raw vector of random number generator state,
##'   returned by `dust_filter_rng_state`
##' @rdname dust_filter_rng_state
##' @export
dust_filter_set_rng_state <- function(filter, rng_state) {
  check_is_dust_filter(filter)
  if (is.null(filter$ptr)) {
    assert_raw_vector(rng_state, length(filter$initial_rng_state))
    filter$initial_rng_state <- rng_state
  } else {
    filter$methods$set_rng_state(filter$ptr, rng_state)
  }
  invisible()
}


check_is_dust_filter <- function(filter, call = parent.frame()) {
  if (!inherits(filter, "dust_filter")) {
    cli::cli_abort("Expected 'filter' to be a 'dust_filter' object",
                   arg = "filter", call = call)
  }
}


## This is something that we should tidy up within mcstate itself,
## there are lots of little utilities we could benefit from.
filter_rng_state <- function(n_particles, n_groups, seed) {
  n_streams <- max(n_groups, 1) * (1 + n_particles)
  mcstate2::mcstate_rng$new(n_streams = n_streams, seed = seed)$state()
}
