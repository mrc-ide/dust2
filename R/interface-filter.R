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
##'   and `time_start` must be no later than than the earliest time.
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
##' @return A `dust_likelihood` object, which can be used with
##'   [dust_likelihood_run]
##'
##' @export
dust_filter_create <- function(generator, time_start, data,
                               n_particles, n_groups = NULL, dt = NULL,
                               ode_control = NULL,
                               n_threads = 1, preserve_group_dimension = FALSE,
                               seed = NULL) {
  call <- environment()
  check_generator_for_filter(generator, "filter", call = call)
  assert_scalar_size(n_particles, allow_zero = FALSE, call = call)
  assert_scalar_logical(preserve_group_dimension, call = call)

  data <- prepare_data(data, n_groups, call = call)
  time_start <- check_time_start(time_start, data$time, call = call)
  time_control <- check_time_control(generator, dt, ode_control, call = call)

  n_groups <- data$n_groups
  preserve_group_dimension <- preserve_group_dimension || n_groups > 1

  n_threads <- check_n_threads(n_threads, n_particles, n_groups)

  if (generator$properties$time_type == "continuous") {
    cli::cli_abort(
      c("Can't use 'dust_filter_create()' with continuous-time models",
        i = paste("We'll support this in future, but you need a place where",
                  "stochasticity enters into the model, and we don't really",
                  "allow for this yet")))
  }

  inputs <- list(time_start = time_start,
                 time = data$time,
                 time_control = time_control,
                 data = data$data,
                 n_particles = n_particles,
                 n_groups = n_groups,
                 n_threads = n_threads,
                 preserve_group_dimension = preserve_group_dimension)

  res <- list2env(
    list(inputs = inputs,
         initialise = filter_create,
         initial_rng_state = filter_rng_state(n_particles, n_groups, seed),
         n_particles = n_particles,
         n_groups = n_groups,
         deterministic = FALSE,
         has_adjoint = FALSE,
         generator = generator,
         methods = generator$methods$filter,
         time_control = time_control,
         preserve_particle_dimension = TRUE,
         preserve_group_dimension = preserve_group_dimension),
    parent = emptyenv())
  class(res) <- c("dust_filter", "dust_likelihood")
  res
}


filter_create <- function(obj, pars) {
  inputs <- obj$inputs
  list2env(
    obj$methods$alloc(pars,
                      inputs$time_start,
                      inputs$time,
                      inputs$time_control,
                      inputs$data,
                      inputs$n_particles,
                      inputs$n_groups,
                      inputs$n_threads,
                      obj$initial_rng_state),
    obj)
  obj$initial_rng_state <- NULL
  obj$packer_state <- monty::monty_packer(array = obj$packing_state)
}


## This is something that we should tidy up within monty itself,
## there are lots of little utilities we could benefit from.
filter_rng_state <- function(n_particles, n_groups, seed) {
  n_streams <- max(n_groups, 1) * (1 + n_particles)
  monty::monty_rng$new(n_streams = n_streams, seed = seed)$state()
}
