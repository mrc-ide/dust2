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
##' @param time A vector of times, each of which has a corresponding
##'   entry in `data`.  The system will stop at each of these times to
##'   compute the likelihood using the compare function.
##'
##' @param data The data to compare against.  This must be a list with
##'   the same length as `time`, each element of which corresponds to
##'   the data required for the system.  If the system is ungrouped then
##'   each element of `data` is a list with elements corresponding to
##'   whatever your system requires.  If your system is grouped, this
##'   should be a list with as many elements as your system has groups,
##'   with each element corresponding to the data your system requires.
##'   We will likely introduce a friendlier data.frame based input
##'   soon.
##'
##' @param n_particles The number of particles to run.  Larger numbers
##'   will give lower variance in the likelihood estimate but run more
##'   slowly.
##'
##' @inheritParams dust_system_create
##' @inheritParams dust_system_simulate
##'
##' @return A `dust_unfilter` object, which can be used with
##'   [dust_unfilter_run]
##'
##' @export
dust_filter_create <- function(generator, time_start, time, data,
                               n_particles, n_groups = 0, dt = 1,
                               index = NULL, seed = NULL) {
  call <- environment()
  check_generator_for_filter(generator, "filter", call = call)
  assert_scalar_size(n_particles, call = call)
  assert_scalar_size(n_groups, allow_zero = TRUE, call = call)
  check_time_sequence(time_start, time, call = call)
  check_dt(dt, call = call)
  check_data(data, length(time), n_groups, call = call)
  check_index(index, call = call)

  create <- function(filter, pars) {
    list2env(
      generator$methods$filter$alloc(pars, time_start, time, dt, data,
                                     n_particles, n_groups, index,
                                     filter$initial_rng_state),
      filter)
    filter$create <- NULL
    filter$initial_rng_state <- NULL
  }
  res <- list2env(
    list(create = create,
         initial_rng_state = filter_rng_state(n_particles, n_groups, seed),
         n_particles = as.integer(n_particles),
         n_groups = as.integer(max(n_groups), 1),
         deterministic = FALSE,
         methods = generator$methods$filter,
         index = index),
    parent = emptyenv())
  class(res) <- "dust_filter"
  res
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
##'   filter was constructed using a non-`NULL` `index` parameter,
##'   the history is restricted to these states.
##'
##' @return A vector of likelihood values, with as many elements as
##'   there are groups.
##'
##' @export
dust_filter_run <- function(filter, pars, initial = NULL,
                            save_history = FALSE) {
  check_is_dust_filter(filter)
  if (is.null(filter$ptr)) {
    if (is.null(pars)) {
      cli::cli_abort("'pars' cannot be NULL, as filter is not initialised",
                     arg = "pars", call = call)
    }
    filter$create(filter, pars)
  } else if (!is.null(pars)) {
    filter$methods$update_pars(filter$ptr, pars, filter$grouped)
  }
  filter$methods$run(filter$ptr, initial, save_history, filter$grouped)
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
  filter$methods$last_history(filter$ptr, filter$grouped)
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


filter_rng_state <- function(n_particles, n_groups, seed) {
  n_streams <- max(n_groups, 1) * (1 + n_particles)
  mcstate2::mcstate_rng$new(n_streams = n_streams, seed = seed)$state()
}
