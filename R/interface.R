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
##' @param pars A list of parameters.  The format of this will depend
##'   on the system.  If `n_groups` is 1 or more, then this must be a
##'   list of length `n_groups` where each element is a list of
##'   parameters for your system.  The default `list()` assumes your
##'   system has no required parameters, but you may need to pass a
##'   list of parameters here.
##'
##' @param n_particles The number of particles to create, defaulting
##'   to a single particle
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
dust_system_create <- function(generator, pars = list(), n_particles = 1,
                               n_groups = 1,
                               time = 0, dt = NULL, ode_control = NULL,
                               seed = NULL, deterministic = FALSE,
                               n_threads = 1,
                               preserve_particle_dimension = FALSE,
                               preserve_group_dimension = FALSE) {
  call <- environment()
  check_is_dust_system_generator(generator, substitute(generator))
  methods <- dust_system_generator_methods(generator)

  time_control <- check_time_control(generator, dt, ode_control, call = call)
  check_time(time, time_control, call = call)

  assert_scalar_size(n_particles, allow_zero = FALSE, call = call)
  assert_scalar_size(n_groups, allow_zero = FALSE, call = call)
  assert_scalar_size(n_threads, allow_zero = FALSE, call = call)
  assert_scalar_logical(preserve_particle_dimension, call = call)
  assert_scalar_logical(preserve_group_dimension, call = call)

  preserve_particle_dimension <- preserve_particle_dimension || n_particles > 1
  preserve_group_dimension <- preserve_group_dimension || n_groups > 1

  pars <- check_pars(pars, n_groups, NULL, preserve_group_dimension)
  res <- methods$alloc(pars, time, time_control, n_particles,
                       n_groups, seed, deterministic, n_threads)

  parameters <- coef(generator)
  if (!is.null(parameters$constant)) {
    parameters <- parameters[!parameters$constant, ]
  }

  ## Here, we augment things slightly
  res$name <- attr(generator, "name")
  res$packer_state <- monty::monty_packer(array = res$packing_state)
  if (!is.null(res$packing_gradient)) {
    res$packer_gradient <-
      monty::monty_packer(array = res$packing_gradient)
  }
  res$n_particles <- as.integer(n_particles)
  res$n_groups <- as.integer(n_groups)
  res$n_threads <- check_n_threads(n_threads, n_particles, n_groups)
  res$deterministic <- deterministic
  res$methods <- methods
  res$properties <- attr(generator, "properties")
  res$parameters <- parameters
  res$preserve_particle_dimension <- preserve_particle_dimension
  res$preserve_group_dimension <- preserve_group_dimension
  res$time_control <- time_control
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
##' length 1 will be recycled as appropriate.  For continuous time
##' systems, we will initialise the solver immediately after setting
##' state, which may cause errors if your initial state is invalid for
##' your system.  There are many ways that you can use this function
##' to set different fractions of state (a subset of states, particles
##' or parameter groups, recycling over any dimensions that are
##' missing).  Please see the Examples section for usage.
##'
##' @title Set system state
##'
##' @inheritParams dust_system_state
##'
##' @param state A matrix or array of state.  If ungrouped, the
##'   dimension order expected is state x particle.  If grouped the
##'   order is state x particle x group.  If you have a grouped system
##'   with 1 particle and `preserve_state_dimension = FALSE` then the
##'   state has size state x group.  You can omit higher dimensions,
##'   so if you pass a vector it will be treated as if all higher
##'   dimensions are length 1 (or if you have a grouped system you
##'   can provide a matrix and treat it as if the third dimension had
##'   length 1).  If you provide any `index_` argument then the length
##'   of the corresponding state dimension must match the index
##'   length.
##'
##' @param index_state An index to control which state variables we
##'   set.  You can use this to set a subset of state variables.
##'
##' @param index_particle An index to control which particles have
##'   their state updated
##'
##' @param index_group An index to control which groups have their
##'   state updated.
##'
##' @return Nothing, called for side effects only
##' @export
##' @examples
##' # Consider a system with 3 particles and 1 group:
##' sir <- dust_example("sir")
##' sys <- dust_system_create(sir(), list(), n_particles = 3)
##' # The state for this system is packed as S, I, R, cases_cumul, cases_inc:
##' dust_unpack_index(sys)
##'
##' # Set all particles to the same state:
##' dust_system_set_state(sys, c(1000, 10, 0, 0, 0))
##' dust_system_state(sys)
##'
##' # We can set everything to different states by passing a vector
##' # with this shape:
##' m <- cbind(c(1000, 10, 0, 0, 0), c(999, 11, 0, 0, 0), c(998, 12, 0, 0, 0))
##' dust_system_set_state(sys, m)
##' dust_system_state(sys)
##'
##' # Or set the state for just one state:
##' dust_system_set_state(sys, 1, index_state = 4)
##' dust_system_state(sys)
##'
##' # If you want to set a different state across particles, you must
##' # provide a *matrix* (a vector always sets the same state into
##' # every particle)
##' dust_system_set_state(sys, rbind(c(1, 2, 3)), index_state = 4)
##' dust_system_state(sys)
##'
##' # This will not work as it can be ambiguous what you are
##' # trying to do:
##' #> dust_system_set_state(sys, c(1, 2, 3), index_state = 4)
##'
##' # State can be set for specific particles:
##' dust_system_set_state(sys, c(900, 100, 0, 0, 0), index_particle = 2)
##' dust_system_state(sys)
##'
##' # And you can combine 'index_particle' with 'index_state' to set
##' # small rectangles of state:
##' dust_system_set_state(sys, matrix(c(1, 2, 3, 4), 2, 2),
##'                       index_particle = 2:3, index_state = 4:5)
##' dust_system_state(sys)
dust_system_set_state <- function(sys, state, index_state = NULL,
                                  index_particle = NULL, index_group = NULL) {
  check_is_dust_system(sys)
  state <- prepare_state(state,
                         index_state,
                         index_particle,
                         index_group,
                         sys$n_state,
                         sys$n_particles,
                         sys$n_groups,
                         sys$preserve_particle_dimension,
                         sys$preserve_group_dimension)
  sys$methods$set_state(sys$ptr, state)
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
  check_time(time, sys$time_control)
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
##'   first time must be no less than the current system time (as
##'   reported by [dust_system_time]).  If your system is discrete,
##'   then times must align to the `dt` used when creating the system.
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
  check_time_sequence(times, sys$time_control)
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
##' @param include_history Boolean, also undocumented.
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
dust_system_internals <- function(sys, include_coefficients = FALSE,
                                  include_history = FALSE) {
  check_is_dust_system(sys)
  if (sys$properties$time_type == "discrete") {
    ## No internals for now, perhaps never?
    return(NULL)
  }
  dat <- sys$methods$internals(sys$ptr, include_coefficients, include_history)
  ret <- data_frame(
    particle = seq_along(dat),
    dydt = I(lapply(dat, "[[", "dydt")),
    step_times = I(lapply(dat, "[[", "step_times")),
    step_size = vnapply(dat, "[[", "step_size"),
    error = vnapply(dat, "[[", "error"),
    n_steps = viapply(dat, "[[", "n_steps"),
    n_steps_accepted = viapply(dat, "[[", "n_steps_accepted"),
    n_steps_rejected = viapply(dat, "[[", "n_steps_rejected"),
    events = I(lapply(dat, function(x) {
      if (is.null(x$events)) NULL else as.data.frame(x$events)
    })))
  if (include_coefficients) {
    ret$coefficients <- I(lapply(dat, "[[", "coefficients"))
  }
  if (include_history) {
    ret$history <- I(lapply(dat, "[[", "history"))
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
  cli::cli_alert_info(format_dimensions(x))
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
  cli::cli_bullets(
    c(i = describe_time(x$properties$time_type, NULL, x$time_control$dt)))
  n_pars <- NROW(x$parameters)
  cli::cli_alert_info(paste(
    "This system has {cli::no(n_pars)} parameter{?s} that can",
    "be updated via {.run dust_system_update_pars}"))
  if (n_pars > 0) {
    cli::cli_bullets(c(">" = "{squote(x$parameters$name)}"))
  }
  cli::cli_alert_info(
    "Use {.help [coef()](stats::coef)} to get more information on parameters")
  invisible(x)
}


##' @export
print.dust_system_generator <- function(x, ...) {
  name <- attr(x, "name")
  properties <- attr(x, "properties")
  default_dt <- attr(x, "default_dt")
  parameters <- attr(x, "parameters")$name
  cli::cli_h1("<dust_system_generator: {name}>")
  if (properties$has_compare) {
    cli::cli_bullets(c(
      i = "This system has 'compare_data' support"))
  }
  cli::cli_bullets(
    c(i = describe_time(properties$time_type, default_dt)))
  n_pars <- length(parameters)
  cli::cli_alert_info("This system has {cli::no(n_pars)} parameter{?s}")
  if (n_pars > 0) {
    cli::cli_bullets(c(">" = "{squote(parameters)}"))
  }
  cli::cli_alert_info(
    paste("Use {.help [dust2::dust_system_create()](dust2::dust_system_create)} to create a system",
          "with this generator"))
  cli::cli_alert_info(
    "Use {.help [coef()](stats::coef)} to get more information on parameters")

  invisible(x)
}


describe_time <- function(time_type, default_dt, dt) {
  if (time_type == "continuous") {
    "This system runs in continuous time"
  } else if (time_type == "discrete") {
    prefix <- "This system runs in discrete time"
    if (rlang::is_missing(dt)) {
      cli::format_inline("{prefix} with a default dt of {default_dt}")
    } else {
      cli::format_inline("{prefix} with dt = {dt}")
    }
  } else if (time_type == "mixed") {
    prefix <- "This system runs in both continuous time and discrete time"
    if (rlang::is_missing(dt)) {
      if (is.null(default_dt)) {
        cli::format_inline("{prefix} with discrete time disabled by default")
      } else {
        cli::format_inline("{prefix} with a default dt of {default_dt}")
      }
    } else {
      if (is.null(dt)) {
        cli::format_inline("{prefix} with discrete time disabled")
      } else {
        cli::format_inline("{prefix} with dt = {dt}")
      }
    }
  }
}

##' @importFrom stats coef
##' @export
coef.dust_system_generator <- function(object, ...) {
  attr(object, "parameters")
}


##' @export
dim.dust_system <- function(x, ...) {
  c(x$n_state,
    if (x$preserve_particle_dimension) x$n_particles,
    if (x$preserve_group_dimension) x$n_groups)
}

##' @export
coef.dust_system <- function(object, ...) {
  object$parameters
}


check_is_dust_system_generator <- function(generator, called_as,
                                           call = parent.frame()) {
  assert_is(generator, "dust_system_generator", call = call)
}


check_is_dust_system <- function(sys, call = parent.frame()) {
  if (!inherits(sys, "dust_system")) {
    cli::cli_abort("Expected 'sys' to be a 'dust_system' object",
                   arg = "sys", call = call)
  }
}


check_time <- function(time, time_control, name = "time",
                       call = parent.frame()) {
  assert_scalar_numeric(time, name = name, call = call)
  dt <- time_control$dt
  if (!is.null(dt) && dt > 0) {
    if (abs(fmod(time, dt)) > sqrt(.Machine$double.eps)) {
      if (dt == 1) {
        cli::cli_abort(
          "'{name}' must be integer-like, because 'dt' is 1",
          arg = name, call = call)
      } else {
        cli::cli_abort(
          "'{name}' must be a multiple of 'dt' ({dt})",
          arg = name, call = call)
      }
    }
  }
}


check_time_sequence <- function(time, time_control,
                                name = deparse(substitute(time)),
                                call = parent.frame()) {
  assert_numeric(time, name = name, call = call)
  if (length(time) == 0) {
    cli::cli_abort("Expected at least one value in '{name}'",
                   arg = name, call = call)
  }

  check_increasing(time, name = name, call = call)

  dt <- time_control$dt
  if ((!is.null(dt)) && (dt <= 1)) {
    #rem <- time %% dt # The problem is 10 %% 0.25 = 0, but 10 %% 0.1 = 0.1
    rem <- fmod(time, dt)
    err <- abs(rem) > sqrt(.Machine$double.eps)
    if (any(err)) {
      i <- which(err)
      detail <- tail_errors(sprintf("'{name}[%d]' (%s) is invalid", i, time[i]))
      if (dt == 1) {
        cli::cli_abort(
          c("Values in '{name}' must be integer-like, because 'dt' is 1",
            detail),
          arg = name, call = call)
      } else {
        cli::cli_abort(
          c("Values in '{name}' must be multiples of 'dt' ({dt})",
            detail),
          arg = name, call = call)
      }
    }
  }

  as.numeric(time)
}


check_increasing <- function(x, name = deparse(substitute(x)),
                             call = parent.frame()) {
  err <- diff(x) <= 0
  if (any(err)) {
    i <- which(err)
    detail <- tail_errors(sprintf(
      "'{name}[%d]' (%s) must be greater than '{name}[%d]' (%s)",
      i + 1, x[i + 1], i, x[i]))
    cli::cli_abort(
      c("Values in '{name}' must be increasing",
        detail),
      arg = name, call = call)
  }
}


format_dimensions <- function(x) {
  if (is.null(x$n_state)) {
    if (x$preserve_group_dimension && x$preserve_particle_dimension) {
      cli::format_inline(paste(
        "{x$n_particles} particle{?s} x {x$n_groups} group{?s}"))
    } else if (x$preserve_group_dimension) {
      cli::format_inline(
        "{x$n_groups} group{?s}")
    } else if (!isFALSE(x$preserve_particle_dimension)) {
      cli::format_inline(
        "{x$n_particles} particle{?s}")
    } else {
      cli::format_inline("single particle")
    }
  } else {
    if (x$preserve_group_dimension && x$preserve_particle_dimension) {
      cli::format_inline(paste(
        "{x$n_state} state x {x$n_particles} particle{?s} x",
        "{x$n_groups} group{?s}"))
    } else if (x$preserve_group_dimension) {
      cli::format_inline(
        "{x$n_state} state x {x$n_groups} group{?s}")
    } else if (!isFALSE(x$preserve_particle_dimension)) {
      cli::format_inline(
        "{x$n_state} state x {x$n_particles} particle{?s}")
    } else {
      cli::format_inline("single particle with {x$n_state} state{?s}")
    }
  }
}


check_time_control <- function(generator, dt, ode_control,
                               call = parent.frame()) {
  dt <- check_system_dt(dt, generator, call = call)
  ode_control <- check_system_ode_control(ode_control, generator, call = call)
  list(dt = dt, ode_control = ode_control)
}


dust_package_env <- function(generator, quiet = FALSE) {
  name <- attr(generator, "name")
  pkg <- attr(generator, "package")
  path <- attr(generator, "path")
  if (is.null(path)) {
    environment(generator)
  } else if (isNamespaceLoaded(pkg)) {
    asNamespace(pkg)
  } else {
    load_temporary_package(path, pkg, quiet)
  }
}


## This does a bunch of bookkeeping to work out if we can set state
## into a system, and works out what we'll need to recycle when
## setting it.
prepare_state <- function(state,
                          index_state,
                          index_particle,
                          index_group,
                          n_state,
                          n_particles,
                          n_groups,
                          preserve_particle_dimension,
                          preserve_group_dimension,
                          name = deparse(substitute(state)),
                          call = parent.frame()) {
  len_from_index <- function(n, idx, name_index = deparse(substitute(idx))) {
    if (is.null(idx)) {
      n
    } else {
      check_index(idx, n, unique = TRUE, name = name_index)
      length(idx)
    }
  }
  len_state <- len_from_index(n_state, index_state)
  len_particles <- len_from_index(n_particles, index_particle)
  len_groups <- len_from_index(n_groups, index_group)

  stopifnot(preserve_particle_dimension || n_particles == 1)
  stopifnot(preserve_group_dimension || n_groups == 1)

  d <- dim2(state)
  rank <- length(d)
  expected <- c(state = len_state,
                particle = if (preserve_particle_dimension) len_particles,
                group = if (preserve_group_dimension) len_groups)
  rank_expected <- length(expected)
  if (rank > rank_expected) {
    cli::cli_abort(
      paste("Expected 'state' to be a {rank_description(rank_expected)}",
            "but was given a {rank_description(rank)}"),
      arg = name, call = call)
  }
  if (rank < rank_expected) {
    d <- c(d, rep(1L, rank_expected - rank))
  }

  ok <- d == expected | c(FALSE, d[-1] == 1)
  if (all(ok)) {
    ## We will access this by position from the C++ code but name it
    ## here for clarity.
    return(list(state = state,
                index_state = index_state,
                index_particle = index_particle,
                index_group = index_group,
                recycle_particle = n_particles > 1 && d[[2]] == 1,
                recycle_group = n_groups > 1 && last(d) == 1))
  }

  if (!ok[[1]]) {
    if (rank == 1) {
      msg <- "Expected '{name}' to have length {len_state}"
    } else if (rank == 2) {
      msg <- "Expected '{name}' to have {len_state} rows"
    } else {
      msg <- "Expected dimension 1 of '{name}' to be length {len_state}"
    }
  } else if (!ok[[2]]) {
    expected_str <-
      if (expected[[2]] == 1) "1" else sprintf("1 or %d", expected[[2]])
    if (rank == 2) {
      msg <- "Expected '{name}' to have {expected_str} columns"
    } else {
      msg <- "Expected dimension 2 of '{name}' to be length {expected_str}"
    }
  } else {
    expected_str <-
      if (expected[[3]] == 1) "1" else sprintf("1 or %d", expected[[3]])
    msg <- "Expected dimension 3 of '{name}' to be length {expected_str}"
  }

  cli::cli_abort(msg, arg = name, call = call)
}


dust_system_generator_methods <- function(generator) {
  name <- attr(generator, "name")
  properties <- attr(generator, "properties")
  env <- dust_package_env(generator)
  time_type <- properties$time_type
  has_compare <- properties$has_compare

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
                    if (time_type != "discrete") "internals")
  methods_compare <- "compare_data"
  methods <- get_methods(c(methods_core, methods_compare), "system", name)
  ok <- !vapply(methods[methods_core], is.null, TRUE)
  stopifnot(all(ok))

  has_compare <- !is.null(methods$compare_data)

  if (has_compare) {
    methods_unfilter <- c("alloc", "run", "update_pars", "last_trajectories",
                          "last_snapshots", "last_state", "last_gradient")
    methods$unfilter <- get_methods(methods_unfilter, "unfilter", name)
    methods_filter <- c("alloc", "run", "update_pars", "last_trajectories",
                        "last_snapshots", "last_state",
                        "rng_state", "set_rng_state")
    methods$filter <- get_methods(methods_filter, "filter", name)
  }

  methods
}


##' @export
"$<-.dust_system" <- function(x, name, value) {
  disable_write("dust_system", name)
}

##' @export
"[[<-.dust_system" <- function(x, i, value) {
  disable_write("dust_system", i)
}

##' @export
"[<-.dust_system" <- function(x, i, ..., value) {
  disable_write("dust_system", i)
}

##' @export
"$<-.dust_system_generator" <- function(x, name, value) {
  disable_write("dust_system_generator", name)
}

##' @export
"[[<-.dust_system_generator" <- function(x, i, value) {
  disable_write("dust_system_generator", i)
}

##' @export
"[<-.dust_system_generator" <- function(x, i, ..., value) {
  disable_write("dust_system_generator", i)
}
