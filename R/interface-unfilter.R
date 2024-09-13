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
##' @return A `dust_likelihood` object, which can be used with
##'   [dust_likelihood_run]
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
         initialise = unfilter_create,
         n_particles = as.integer(n_particles),
         n_groups = as.integer(n_groups),
         deterministic = TRUE,
         has_adjoint = generator$properties$has_adjoint,
         generator = generator,
         methods = generator$methods$unfilter,
         index_state = index_state,
         preserve_particle_dimension = preserve_particle_dimension,
         preserve_group_dimension = preserve_group_dimension),
    parent = emptyenv())
  class(res) <- c("dust_unfilter", "dust_likelihood")
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
