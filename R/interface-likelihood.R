##' Compute a log likelihood based on a dynamical model, either from
##' particle filter (created with [dust_filter_create]) or a
##' deterministic model (created with [dust_unfilter_create]).
##'
##' @title Compute likelihood
##'
##' @param obj A `dust_filter` object, created by
##'   [dust_filter_create] or a `dust_unfilter` object created by
##'   [dust_unfilter_create]
##'
##' @param pars Optional parameters to compute the likelihood with.
##'   If not provided, parameters are not updated
##'
##' @param initial Optional initial conditions, as a matrix (state x
##'   particle) or 3d array (state x particle x group).  If not
##'   provided, the system initial conditions are used.
##'
##' @param save_trajectories Logical, indicating if the trajectories
##'   from the simulation history should be saved while the simulation
##'   runs; this has a small overhead in runtime and in memory.
##'   Particle trajectories will be saved at each time for which you
##'   have data.  If `index_state` is non-`NULL`, the trajectories are
##'   limited to these states.
##'
##' @param adjoint Optional logical, indicating if we should enable
##'   adjoint history saving. This is enabled by default if your model
##'   has an adjoint, but can be disabled or enabled even when your
##'   model does not support adjoints! But if you don't actually have
##'   an adjoint you will not be able to compute gradients.  This has
##'   no effect for stochastic models.
##'
##' @param index_state An optional vector of state indices to save
##'   where `save_trajectories` is `TRUE`.  If omitted all state is saved
##'   and if `save_trajectories = FALSE` this has no effect.
##'
##' @param index_group An optional vector of group indices to run the
##'   calculation for.  You can use this to run a subset of possible
##'   groups, once `obj` is initialised (this argument must be `NULL`
##'   on the **first** call).
##'
##' @return A vector of likelihood values, with as many elements as
##'   there are groups.
##'
##' @export
dust_likelihood_run <- function(obj, pars, initial = NULL,
                                save_trajectories = FALSE,
                                adjoint = NULL,
                                index_state = NULL,
                                index_group = NULL) {
  check_is_dust_likelihood(obj)
  index_group <- check_index(index_group, max = obj$n_groups,
                             unique = TRUE)
  if (!is.null(pars)) {
    pars <- check_pars(pars, obj$n_groups, index_group,
                       obj$preserve_group_dimension)
  }

  if (is.null(obj$ptr)) {
    if (is.null(pars)) {
      cli::cli_abort("'pars' cannot be NULL, as 'obj' is not initialised",
                     arg = "pars")
    }
    if (!is.null(index_group)) {
      cli::cli_abort(
        "'index_group' must be NULL, as 'obj' is not initialised",
        arg = "index_group")
    }
    obj$initialise(obj, pars)
  } else if (!is.null(pars)) {
    obj$methods$update_pars(obj$ptr, pars, index_group)
  }
  if (is.null(adjoint)) {
    adjoint <- obj$has_adjoint
  } else {
    assert_scalar_logical(adjoint, call = environment())
    ## Here we might check for adjoint = TRUE in a stochastic model?
  }

  ## This must be evaluated *after* we're guaranteed to be initialised
  ## so that we can access `n_groups`
  index_state <- check_index(index_state, max = obj$n_state,
                             unique = TRUE)
  if (!is.null(initial)) {
    initial <- prepare_state(initial,
                             NULL, # index_state
                             NULL, # index_particle
                             index_group,
                             obj$n_state,
                             obj$n_particles,
                             obj$n_groups,
                             obj$preserve_particle_dimension,
                             obj$preserve_group_dimension)
  }
  obj$methods$run(obj$ptr,
                  initial,
                  save_trajectories,
                  adjoint,
                  index_state,
                  index_group,
                  obj$preserve_particle_dimension,
                  obj$preserve_group_dimension)
}


dust_likelihood_ensure_initialised <- function(obj, pars) {
  is_uninitialised <- is.null(obj$ptr)
  if (is_uninitialised) {
    index_group <- NULL
    pars <- check_pars(pars, obj$n_groups, index_group,
                       obj$preserve_group_dimension)
    obj$initialise(obj, pars)
  }
  is_uninitialised
}


##' Create an independent copy of a likelihood object.  The new object
##' is decoupled from the random number streams of the parent object.
##' It is also decoupled from the *state size* of the parent object,
##' so you can use this to create a new object where the system is
##' fundamentally different but everything else is the same.
##'
##' @title Create copy of a dust likelihood object
##'
##' @inheritParams dust_likelihood_run
##'
##' @param seed The seed for the particle filter (see
##'   [dust_filter_create])
##'
##' @return A new `dust_likelihood` object
dust_likelihood_copy <- function(obj, seed = NULL) {
  check_is_dust_likelihood(obj)
  dst <- new.env(parent = emptyenv())
  for (nm in setdiff(names(obj), "initial_rng_state")) {
    dst[[nm]] <- obj[[nm]]
  }
  if (!is.null(obj$initial_rng_state)) {
    dst$initial_rng_state <-
      filter_rng_state(obj$n_particles, obj$n_groups, seed)
  }
  class(dst) <- class(obj)
  dst
}


##' Fetch the last trajectories created by running a likelihood.  This
##' errors if the last call to [dust_likelihood_run] did not use
##' `save_trajectories = TRUE`.  We return the states and groups that
##' were run via the `index_state` and `index_group` arguments to
##' [dust_likelihood_run].
##'
##' @title Fetch last likelihood trajectories
##'
##' @inheritParams dust_likelihood_run
##'
##' @param select_random_particle Logical, indicating if we should
##'   return a history for one randomly selected particle (rather than
##'   the entire set of trajectories over all particles).  If this is
##'   `TRUE`, the particle will be selected independently for each
##'   group, if the object is grouped.  This option is intended to
##'   help select a representative trajectory during an MCMC.  When
##'   `TRUE`, we drop the `particle` dimension of the return value.
##'
##' @return An array.  If ungrouped this will have dimensions `state`
##'   x `particle` x `time`, and if grouped then `state` x `particle`
##'   x `group` x `time`.  If `select_random_particle = TRUE`, the
##'   second (particle) dimension will be dropped.
##'
##' @export
dust_likelihood_last_trajectories <- function(obj,
                                              select_random_particle = FALSE) {
  check_is_dust_likelihood(obj)
  if (is.null(obj$ptr)) {
    cli::cli_abort(c(
      "Trajectories are not current",
      i = "Likelihood has not yet been run"))
  }
  assert_scalar_logical(select_random_particle)
  obj$methods$last_trajectories(obj$ptr,
                                select_random_particle,
                                obj$preserve_particle_dimension,
                                obj$preserve_group_dimension)
}


##' Get the last state from a likelihood.
##'
##' @title Get likelihood state
##'
##' @inheritParams dust_likelihood_last_trajectories
##'
##' @return An array.  If ungrouped this will have dimensions `state`
##'   x `particle`, and if grouped then `state` x `particle` x
##'   `group`.  If `select_random_particle = TRUE`, the second
##'   (particle) dimension will be dropped.  This is the same as the
##'   state returned by [dust_likelihood_last_trajectories] without the time
##'   dimension but also without any state index applied (i.e., we
##'   always return all state).
##'
##' @export
dust_likelihood_last_state <- function(obj, select_random_particle = FALSE) {
  check_is_dust_likelihood(obj)
  if (is.null(obj$ptr)) {
    cli::cli_abort(c(
      "State is not current",
      i = "Likelihood has not yet been run"))
  }
  assert_scalar_logical(select_random_particle)
  obj$methods$last_state(obj$ptr,
                         select_random_particle,
                         obj$preserve_particle_dimension,
                         obj$preserve_group_dimension)
}


##' Get random number generator (RNG) state from the particle filter.
##'
##' @title Get filter RNG state
##'
##' @inheritParams dust_likelihood_run
##'
##' @return A raw vector, this could be quite long.  Later we will
##'   describe how you might reseed a filter or system with this state.
##'
##' @export
dust_likelihood_rng_state <- function(obj) {
  check_is_dust_likelihood(obj)
  if (is.null(obj$ptr)) {
    obj$initial_rng_state
  } else {
    obj$methods$rng_state(obj$ptr)
  }
}


##' @param rng_state A raw vector of random number generator state,
##'   returned by `dust_likelihood_rng_state`
##'
##' @rdname dust_likelihood_rng_state
##'
##' @export
dust_likelihood_set_rng_state <- function(obj, rng_state) {
  check_is_dust_likelihood(obj)
  if (check_externalptr(obj$ptr, TRUE)) {
    obj$methods$set_rng_state(obj$ptr, rng_state)
  } else {
    assert_raw(rng_state, length(obj$initial_rng_state))
    assign("initial_rng_state", rng_state, obj)
  }
  invisible()
}


##' Fetch the last gradient created by running an likelihood.  This
##' errors if the last call to [dust_likelihood_run] did not use
##' `adjoint = TRUE`.  The first time you call this (after a
##' particular set of parameters) it will trigger running the reverse
##' model.
##'
##' @title Fetch last likelihood gradient
##'
##' @inheritParams dust_likelihood_run
##'
##' @return A vector (if ungrouped) or a matrix (if grouped).
##'
##' @export
dust_likelihood_last_gradient <- function(obj) {
  check_is_dust_likelihood(obj)
  if (is.null(obj$ptr)) {
    cli::cli_abort(c(
      "Gradient is not current",
      i = "Likelihood has not yet been run"))
  }
  obj$methods$last_gradient(obj$ptr,
                            obj$preserve_particle_dimension,
                            obj$preserve_group_dimension)
}


check_is_dust_likelihood <- function(obj, name = deparse(substitute(obj)),
                                     arg = name, call = parent.frame()) {
  assert_is(obj, "dust_likelihood", name = name, arg = arg, call = call)
}


##' @export
print.dust_likelihood <- function(x, ...) {
  cli::cli_h1("<dust_likelihood ({attr(x$generator, 'name')})>")
  cli::cli_alert_info(format_dimensions(x))
  type <- if (x$deterministic) "deterministic" else "stochastic"
  cli::cli_alert_info("The likelihood is {type}")
  ## TODO:
  ## * has gradient
  ## Link to docs
  time_type <- attr(x$generator, "properties")$time_type
  cli::cli_bullets(
    c(i = describe_time(time_type, NULL, x$time_control$dt)))
  invisible(x)
}
