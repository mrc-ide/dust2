##' Create a system generator.  This function is not intended to be
##' used by users directly, but will be called from packages built
##' using [dust_package]
##'
##' @title Create a system generator
##'
##' @param name The name of the generator
##'
##' @param time_type The time type (discrete or continuous).  Using
##'   the wrong time here will lead to crashes or failure to create
##'   the generator.
##'
##' @param env The environment where the generator is defined.
##'
##' @return A `dust_system_generator` object
##'
##' @export
##' @keywords internal
##' @examples
##' # This creates the "sir" generator
##' dust_system_generator("sir", "discrete", asNamespace("dust2"))
##'
##' # This is the same code as in "dust2:::sir", except there we find
##' # the correct environment automatically
##' dust2:::sir
dust_system_generator <- function(name, time_type,
                                  env = parent.env(parent.frame())) {
  ## I don't love that this requires running through sprintf() each
  ## time we create a generator, but using a function for the generator (see
  ## sir()), rather than an object, means that it's easier to think
  ## about the dependencies among packages.  This is also essentially
  ## how DBI works.
  get_methods <- function(nms, component, name) {
    set_names(
      lapply(sprintf("dust2_%s_%s_%s", component, name, nms),
             function(x) env[[x]]),
      nms)
  }

  methods_core <- c("alloc",
                    "state", "set_state", "set_state_initial",
                    "time", "set_time",
                    "rng_state", "set_rng_state",
                    "update_pars",
                    "run_to_time", "simulate",
                    "reorder",
                    if (time_type == "continuous") "internals")
  methods_compare <- "compare_data"
  methods <- get_methods(c(methods_core, methods_compare), "system", name)
  ok <- !vapply(methods[methods_core], is.null, TRUE)
  stopifnot(all(ok))

  has_compare <- !is.null(methods$compare_data)

  if (has_compare) {
    methods_unfilter <- c("alloc", "run", "update_pars", "last_history",
                          "last_gradient")
    methods$unfilter <- get_methods(methods_unfilter, "unfilter", name)
    methods_filter <- c("alloc", "run", "update_pars", "last_history",
                        "rng_state", "set_rng_state")
    methods$filter <- get_methods(methods_filter, "filter", name)
  }

  properties <- list(
    time_type = time_type,
    has_compare = has_compare,
    has_adjoint = !is.null(methods$unfilter$last_gradient))

  ret <- list(name = name,
              methods = methods,
              properties = properties)

  class(ret) <- "dust_system_generator"
  ret
}


##' Create a dust system object from a system generator.  This allocates a
##' system and sets an initial set of parameters.  Once created you can use
##' other dust functions to interact with it.
##'
##' # Parallelisation
##'
##' Many calculations within a dust system can be parallelised
##' straightforwardly - the most important of these is typically
##' running the model (via [dust_system_run_to_time] or
##' [dust_system_simulate]) but we also parallelise
##' [dust_system_set_state_initial], [dust_system_compare_data] and
##' even [dust_system_reorder].  You need to set the number of threads
##' for parallelism at system creation, and this number cannot be
##' usefully larger than `n_particles` (or `n_particles * n_groups` if
##' you have a grouped system).
##'
##' @title Create a dust system object
##'
##' @param generator A system generator object, with class
##'   `dust_system_generator`
##'
##' @param pars A list of parameters.  The format of this will depend on the
##'   system.  If `n_groups` is 1 or more, then this must be a list of length
##'   `n_groups` where each element is a list of parameters for your system.
##'
##' @param n_particles The number of particles to create.
##'
##' @param n_groups Optionally, the number of parameter groups
##'
##' @param time The initial time, defaults to 0
##'
##' @param dt The time step for discrete time systems, defaults to 1
##'   if not given.  It is an error to provide a non-NULL argument
##'   with continuous-time systems.
##'
##' @param ode_control The ODE integration control for continuous time
##'   systems.  Defaults to the default return of [dust_ode_control].
##'   It is an error to provide this with discrete-time systems.
##'
##' @param seed Optionally, a seed.  Otherwise we respond to R's RNG seed on
##'   initialisation.
##'
##' @param deterministic Logical, indicating if the system should be
##'   allocated in deterministic mode.
##'
##' @param n_threads Integer, the number of threads to use in
##'   parallelisable calculations.  See Details.
##'
##' @param preserve_particle_dimension Logical, indicating if output
##'   from the system should preserve the particle dimension in the
##'   case where a single particle is run.  In the case where more
##'   than one particle is run, this argument has no effect as the
##'   dimension is always preserved.
##'
##' @param preserve_group_dimension Logical, indicating if state and
##'   output from the system should preserve the group dimension in
##'   the case where a single group is run.  In the case where more
##'   than one group is run, this argument has no effect as the
##'   dimension is always preserved.
##'
##' @return A `dust_system` object, with opaque format.
##'
##' @export
dust_system_create <- function(generator, pars, n_particles, n_groups = 1,
                               time = 0, dt = NULL, ode_control = NULL,
                               seed = NULL, deterministic = FALSE,
                               n_threads = 1,
                               preserve_particle_dimension = FALSE,
                               preserve_group_dimension = FALSE) {
  call <- environment()
  check_is_dust_system_generator(generator, substitute(generator))
  ## check_time(time, call = call)
  assert_scalar_size(n_particles, allow_zero = FALSE, call = call)
  assert_scalar_size(n_groups, allow_zero = FALSE, call = call)
  assert_scalar_size(n_threads, allow_zero = FALSE, call = call)
  assert_scalar_logical(preserve_particle_dimension, call = call)
  assert_scalar_logical(preserve_group_dimension, call = call)

  preserve_particle_dimension <- preserve_particle_dimension || n_particles > 1
  preserve_group_dimension <- preserve_group_dimension || n_groups > 1

  pars <- check_pars(pars, n_groups, NULL, preserve_group_dimension)
  if (generator$properties$time_type == "discrete") {
    if (!is.null(ode_control)) {
      cli::cli_abort("Can't use 'ode_control' with discrete-time systems")
    }
    dt <- check_dt(dt %||% 1, call = call)
    res <- generator$methods$alloc(pars, time, dt, n_particles,
                                   n_groups, seed, deterministic,
                                   n_threads)
  } else {
    if (!is.null(dt)) {
      cli::cli_abort("Can't use 'dt' with continuous-time systems")
    }
    if (is.null(ode_control)) {
      ode_control <- dust_ode_control()
    } else {
      assert_is(ode_control, "dust_ode_control", call = environment())
    }
    res <- generator$methods$alloc(pars, time, ode_control, n_particles,
                                   n_groups, seed, deterministic, n_threads)
  }
  ## Here, we augment things slightly
  res$name <- generator$name
  res$n_particles <- as.integer(n_particles)
  res$n_groups <- as.integer(n_groups)
  res$n_threads <- check_n_threads(n_threads, n_particles, n_groups)
  res$deterministic <- deterministic
  res$methods <- generator$methods
  res$properties <- generator$properties
  res$preserve_particle_dimension <- preserve_particle_dimension
  res$preserve_group_dimension <- preserve_group_dimension
  res$ode_control <- ode_control
  res$dt <- dt
  class(res) <- "dust_system"
  res
}


##' Extract system state
##'
##' @title Extract system state
##'
##' @param sys A `dust_system` object
##'
##' @param index_state Index of the state to fetch, if you would like
##'   only a subset
##'
##' @param index_particle Index of the particle to fetch, if you would
##'   like a subset
##'
##' @param index_group Index of the group to fetch, if you would like
##'   a subset
##'
##' @return An array of system state.  If your system is ungrouped
##'   (i.e., `n_groups = 1` and `preserve_group_dimension = FALSE`),
##'   then this has two dimensions (state, particle).  If grouped,
##'   this has three dimensions (state, particle, group)
##'
##' @seealso [dust_system_set_state()] for setting state and
##'   [dust_system_set_state_initial()] for setting state to the
##'   system-specific initial conditions.
##'
##' @export
dust_system_state <- function(sys, index_state = NULL, index_particle = NULL,
                              index_group = NULL) {
  check_is_dust_system(sys)
  sys$methods$state(sys$ptr, index_state, index_particle, index_group,
                    sys$preserve_particle_dimension,
                    sys$preserve_group_dimension)
}


##' Set system state.  Takes a multidimensional array (2- or 3d
##' depending on if the system is grouped or not).  Dimensions of
##' length 1 will be recycled as appropriate.
##'
##' @title Set system state
##'
##' @inheritParams dust_system_state
##'
##' @param state A matrix or array of state.  If ungrouped, the
##'   dimension order expected is state x particle.  If grouped the
##'   order is state x particle x group.
##'
##' @return Nothing, called for side effects only
##' @export
dust_system_set_state <- function(sys, state) {
  check_is_dust_system(sys)
  ## TODO: check rank etc here (mrc-5565), and support
  ## preserve_particle_dimension
  sys$methods$set_state(sys$ptr, state, sys$preserve_group_dimension)
  invisible()
}


##' Set system state from a system's initial conditions.  This may depend
##' on the current time.
##'
##' @title Set system state to initial conditions
##'
##' @inheritParams dust_system_state
##'
##' @return Nothing, called for side effects only
##' @export
dust_system_set_state_initial <- function(sys) {
  check_is_dust_system(sys)
  sys$methods$set_state_initial(sys$ptr)
  invisible()
}


##' Fetch the current time from the system.
##'
##' @title Fetch system time
##'
##' @inheritParams dust_system_state
##'
##' @return A single numeric value
##' @seealso [dust_system_set_time]
##' @export
dust_system_time <- function(sys) {
  check_is_dust_system(sys)
  sys$methods$time(sys$ptr)
}


##' Set time into the system.  This updates the time to the provided
##' value but does not affect the state.  You may want to call
##' [dust_system_set_state] or [dust_system_set_state_initial] after
##' calling this.
##'
##' @title Set system time
##'
##' @inheritParams dust_system_time
##'
##' @param time The time to set.  Currently this must be an
##'   integer-like value, but in future we will allow setting to any
##'   multiple of `dt`.
##'
##' @return Nothing, called for side effects only
##' @export
dust_system_set_time <- function(sys, time) {
  check_is_dust_system(sys)
  ## check_time(time) # TODO
  sys$methods$set_time(sys$ptr, time)
  invisible()
}


##' Fetch, and set, the random number generator (RNG) state from the
##' system.
##'
##' @title Fetch and set rng state
##'
##' @inheritParams dust_system_state
##'
##' @return A raw vector, this could be quite long.
##'
##' @seealso You can pass the state you get back from this function as
##'   the seed object to `dust_system_create` and
##'   `dust_system_set_rng_state`
##'
##' @export
dust_system_rng_state <- function(sys) {
  check_is_dust_system(sys)
  sys$methods$rng_state(sys$ptr)
}


##' @param rng_state A raw vector of random number generator state,
##'   returned by `dust_system_rng_state`
##' @rdname dust_system_rng_state
##' @export
dust_system_set_rng_state <- function(sys, rng_state) {
  check_is_dust_system(sys)
  sys$methods$set_rng_state(sys$ptr, rng_state)
  invisible()
}


##' Update parameters used by the system.  This can be used to update a
##' subset of parameters that do not change the extent of the system,
##' and will be potentially faster than creating a new system object.
##'
##' @title Update parameters
##'
##' @inheritParams dust_system_state
##'
##' @param pars Parameters to set into the system.
##'
##' @return Nothing, called for side effects only
##' @export
dust_system_update_pars <- function(sys, pars) {
  check_is_dust_system(sys)
  pars <- check_pars(pars, sys$n_groups, NULL, sys$preserve_group_dimension)
  sys$methods$update_pars(sys$ptr, pars)
  invisible()
}


##' Run a system, advancing time and the state by repeatedly running its
##' `update` method.  You can advance a system up to a time (which must be in
##' the future).
##'
##' @title Run system
##'
##' @inheritParams dust_system_state
##'
##' @param time Time to run to
##'
##' @return Nothing, called for side effects only
##'
##' @export
dust_system_run_to_time <- function(sys, time) {
  check_is_dust_system(sys)
  sys$methods$run_to_time(sys$ptr, time)
  invisible()
}


##' Simulate a system over a series of times, returning an array of
##' output.  This output can be quite large, so you may filter states
##' according to some index.
##'
##' @title Simulate system
##'
##' @inheritParams dust_system_state
##'
##' @param times A vector of times.  They must be increasing, and the
##'   first time must be no less than the current system time
##'   (as reported by [dust_system_time])
##'
##' @param index_state An optional index of states to extract.  If
##'   given, then we subset the system state on return.  You can use
##'   this to return fewer system states than the system ran with, to
##'   reorder states, or to name them on exit (names present on the
##'   index will be copied into the rownames of the returned array).
##'
##' @return An array with 3 dimensions (state x particle x time) or 4
##'   dimensions (state x particle x group x time) for a grouped
##'   system.
##'
##' @export
dust_system_simulate <- function(sys, times, index_state = NULL) {
  check_is_dust_system(sys)
  ## TODO: check time sequence, except for the first?
  ret <- sys$methods$simulate(sys$ptr,
                              times,
                              index_state,
                              sys$preserve_particle_dimension,
                              sys$preserve_group_dimension)
  if (!is.null(index_state) && !is.null(names(index_state))) {
    rownames(ret) <- names(index_state)
  }
  ret
}


##' Reorder states within a system.  This function is primarily used
##' for debugging and may be removed from the interface if it is not
##' generally useful.
##'
##' @title Reorder states
##'
##' @inheritParams dust_system_state
##'
##' @param index The parameter ordering.  For an ungrouped system this
##'   is a vector where each element is the parameter index (if
##'   element `i` is `j` then after reordering the `i`th particle will
##'   have the state previously used by `j`).  All elements must lie
##'   in `[1, n_particles]`, repetition of an index is allowed (so
##'   that many new particles may have the state as one old particle).
##'   If the system is grouped, `index` must be a matrix with
##'   `n_particles` rows and `n_groups` columns, with each column
##'   corresponding to the reordering for a group.
##'
##' @return Nothing, called for side effects only.
##'
##' @export
dust_system_reorder <- function(sys, index) {
  check_is_dust_system(sys)
  ## TODO: check reorder dimensionality? less important here because
  ## this is a debugging method.
  sys$methods$reorder(sys$ptr, index)
  invisible()
}


##' Return internal data from the system.  This is intended for
##' debugging only, and all formats are subject to change.
##'
##' @title Fetch system internals
##'
##' @inheritParams dust_system_state
##'
##' @param include_coefficients Boolean, indicating if interpolation
##'   coefficients should be included in the output.  These are
##'   intentionally undocumented for now.
##'
##' @return If `sys` is a discrete-time system, this function returns
##'   `NULL`, as no internal data is stored.  Otherwise, for a
##'   continuous-time system we return a `data.frame` of statistics
##'   with one row per particle.  Most of the columns are simple
##'   integers or numeric values, but `dydt` (the current derivative
##'   of the target function with respect to time) and `step_times`
##'   (times that the solver has stopped at, if
##'   `debug_record_step_times` is in [dust_ode_control] was set to
##'   `TRUE`) will be a list of columns, each element of which is a
##'   numeric vector.  If `include_coefficients` is `TRUE`, the
##'   `coefficients` column exists and holds a list of coefficients
##'   (the structure of these may change over time, too).
##'
##' @export
dust_system_internals <- function(sys, include_coefficients = FALSE) {
  check_is_dust_system(sys)
  if (sys$properties$time_type == "discrete") {
    ## No internals for now, perhaps never?
    return(NULL)
  }
  dat <- sys$methods$internals(sys$ptr, include_coefficients)
  ret <- data_frame(
    particle = seq_along(dat),
    dydt = I(lapply(dat, "[[", "dydt")),
    step_times = I(lapply(dat, "[[", "step_times")),
    step_size = vnapply(dat, "[[", "step_size"),
    error = vnapply(dat, "[[", "error"),
    n_steps = viapply(dat, "[[", "n_steps"),
    n_steps_accepted = viapply(dat, "[[", "n_steps_accepted"),
    n_steps_rejected = viapply(dat, "[[", "n_steps_rejected"))
  if (include_coefficients) {
    ret$coefficients <- I(lapply(dat, "[[", "coefficients"))
  }
  ret
}


##' Compare current system state against data.  This is only supported
##' for systems that have 'compare_data' support (i.e., the system
##' definition includes a `compare_data` method).  The current state
##' in the system ([dust_system_state]) is compared against the data
##' provided as `data`.
##'
##' @title Compare system state against data
##'
##' @inheritParams dust_system_state
##'
##' @param data The data to compare against. If the system is ungrouped
##'   then `data` is a list with elements corresponding to whatever
##'   your system requires.  If your system is grouped, this should be a
##'   list with as many elements as your system has groups, with each
##'   element corresponding to the data your system requires.
##'
##' @return A numeric vector with as many elements as your system has
##'   groups, corresponding to the log-likelihood of the data for each
##'   group.
##'
##' @export
dust_system_compare_data <- function(sys, data) {
  check_is_dust_system(sys)
  if (!sys$properties$has_compare) {
    ## This moves into something general soon?
    cli::cli_abort(
      paste("Can't compare against data; the '{sys$name}' system does not",
            "have 'compare_data' support"),
      arg = "system")
  }
  if (!sys$preserve_group_dimension) {
    data <- list(data)
  }
  sys$methods$compare_data(sys$ptr, data, sys$preserve_particle_dimension,
                           sys$preserve_group_dimension)
}


##' @export
print.dust_system <- function(x, ...) {
  cli::cli_h1("<dust_system: {x$name}>")
  ## TODO: Special treatment if not preserve particle dimension
  if (x$preserve_group_dimension) {
    cli::cli_alert_info(paste(
      "{x$n_state} state x {x$n_particles} particle{?s} x",
      "{x$n_groups} group{?s}"))
  } else {
    cli::cli_alert_info("{x$n_state} state x {x$n_particles} particle{?s}")
  }
  if (x$deterministic) {
    cli::cli_bullets(c(
      i = "This system is deterministic"))
  }
  if (x$properties$has_compare) {
    cli::cli_bullets(c(
      i = "This system has 'compare_data' support"))
  }
  if (x$properties$has_adjoint) {
    cli::cli_bullets(c(
      i = "This system has 'adjoint' support, and can compute gradients"))
  }
  cli::cli_bullets(c(
    i = "This system runs in {x$properties$time_type} time"))
  ## Later, we might print some additional capabilities of the system
  ## here, such as if it can be used with a filter, a summary of its
  ## parameters (once we know how to access that), etc.
  invisible(x)
}


##' @export
print.dust_system_generator <- function(x, ...) {
  cli::cli_h1("<dust_system_generator: {x$name}>")
  ## Later, we might print some additional capabilities of the system
  ## here, such as if it can be used with a filter, a summary of its
  ## parameters (once we know how to access that), etc.
  cli::cli_alert_info(
    "Use 'dust2::dust_system_create()' to create a system with this generator")
  if (x$properties$has_compare) {
    cli::cli_bullets(c(
      i = "This system has 'compare_data' support"))
  }
  cli::cli_bullets(c(
    i = "This system runs in {x$properties$time_type} time"))
  invisible(x)
}


##' @export
dim.dust_system <- function(x, ...) {
  c(x$n_state,
    if (x$preserve_particle_dimension) x$n_particles,
    if (x$preserve_group_dimension) x$n_groups)
}


check_is_dust_system_generator <- function(generator, called_as,
                                           call = parent.frame()) {
  if (!inherits(generator, "dust_system_generator")) {
    hint <- NULL
    if (is_uncalled_generator(generator) && is.symbol(called_as)) {
      hint <- c(
        i = "Did you mean '{deparse(called_as)}()' (i.e., with parentheses)")
    }
    cli::cli_abort(
      c("Expected 'generator' to be a 'dust_system_generator' object",
        hint),
      arg = "generator")
  }
}


check_is_dust_system <- function(sys, call = parent.frame()) {
  if (!inherits(sys, "dust_system")) {
    cli::cli_abort("Expected 'sys' to be a 'dust_system' object",
                   arg = "sys", call = call)
  }
}


is_uncalled_generator <- function(sys) {
  if (!is.function(sys)) {
    return(FALSE)
  }
  code <- body(sys)
  rlang::is_call(code, "{") &&
    length(code) == 2 &&
    rlang::is_call(code[[2]], "dust_system")
}
