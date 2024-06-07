dust_system_generator <- function(name, env = parent.env(parent.frame())) {
  prefix <- sprintf("dust2_discrete_%s", name)
  ## I don't love that this requires running through sprintf() each
  ## time we create a generator, but using a function for the generator (see
  ## sir()), rather than an object, means that it's easier to think
  ## about the dependencies among packages.  This is also essentially
  ## how DBI works.
  get_methods <- function(nms, prefix) {
    set_names(
      lapply(sprintf("dust2_discrete_%s_%s", prefix, nms),
             function(x) env[[x]]),
      nms)
  }

  methods_core <- c("alloc",
                    "state", "set_state", "set_state_initial",
                    "time", "set_time",
                    "rng_state", "set_rng_state",
                    "update_pars",
                    "run_steps", "run_to_time", "simulate",
                    "reorder")
  methods_compare <- "compare_data"
  methods <- get_methods(c(methods_core, methods_compare), name)
  ok <- !vapply(methods[methods_core], is.null, TRUE)
  stopifnot(all(ok))

  has_compare <- !is.null(methods$compare_data)

  if (has_compare) {
    methods_unfilter <- c("alloc", "run", "update_pars", "last_history",
                          "last_gradient")
    methods$unfilter <-
      get_methods(methods_unfilter, sprintf("%s_unfilter", name))
    methods_filter <- c("alloc", "run", "update_pars", "last_history",
                        "rng_state", "set_rng_state")
    methods$filter <-
      get_methods(methods_filter, sprintf("%s_filter", name))
  }

  properties <- list(
    has_compare = has_compare,
    has_adjoint = !is.null(methods$unfilter$last_gradient))

  ret <- list(name = name,
              methods = methods,
              properties = properties)
  ## TODO: check that alloc exists, then go through and add
  ## properties.
  class(ret) <- "dust_system_generator"
  ret
}


##' Create a dust system object from a system generator.  This allocates a
##' system and sets an initial set of parameters.  Once created you can use
##' other dust functions to interact with it.
##'
##' @title Create a dust object
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
##' @param dt The time step for the system, defaults to 1
##'
##' @param seed Optionally, a seed.  Otherwise we respond to R's RNG seed on
##'   initialisation.
##'
##' @param deterministic Logical, indicating if the system should be
##'   allocated in deterministic mode.
##'
##' @return A `dust_system` object, with opaque format.
##'
##' @export
dust_system_create <- function(generator, pars, n_particles, n_groups = 0,
                               time = 0, dt = 1,
                               seed = NULL, deterministic = FALSE) {
  check_is_dust_system_generator(generator, substitute(generator))
  res <- generator$methods$alloc(pars, time, dt, n_particles, n_groups,
                                 seed, deterministic)
  ## Here, we augment things slightly
  res$name <- generator$name
  res$n_particles <- as.integer(n_particles)
  res$n_groups <- as.integer(max(n_groups), 1)
  res$deterministic <- deterministic
  res$methods <- generator$methods
  res$properties <- generator$properties
  class(res) <- "dust_system"
  res
}


##' Extract system state
##'
##' @title Extract system state
##'
##' @param sys A `dust_system` object
##'
##' @return An array of system state.  If your system is ungrouped, then
##'   this has two dimensions (state, particle).  If grouped, this has
##'   three dimensions (state, particle, group)
##'
##' @seealso [dust_system_set_state()] for setting state and
##'   [dust_system_set_state_initial()] for setting state to the
##'   system-specific initial conditions.
##'
##' @export
dust_system_state <- function(sys) {
  check_is_dust_system(sys)
  sys$methods$state(sys$ptr, sys$grouped)
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
  sys$methods$set_state(sys$ptr, state, sys$grouped)
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
  sys$methods$update_pars(sys$ptr, pars, sys$grouped)
  invisible()
}


##' Run a system, advancing time and the state by repeatedly running
##' its `update` method.  You can advance a system either a fixed
##' (positive) number of steps, or up to a time (which must be in the
##' future).
##'
##' @title Run system
##'
##' @inheritParams dust_system_state
##'
##' @param steps The number of steps to run forward
##'
##' @return Nothing, called for side effects only
##'
##' @export
##'
##' @rdname dust_system_run
dust_system_run_steps <- function(sys, steps) {
  check_is_dust_system(sys)
  sys$methods$run_steps(sys$ptr, steps)
  invisible()
}


##' @param time Time to run to
##' @export
##' @rdname dust_system_run
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
##' @param index An optional index of states to extract.  If given,
##'   then we subset the system state on return.  You can use this to
##'   return fewer system states than the system ran with, to reorder
##'   states, or to name them on exit (names present on the index will
##'   be copied into the rownames of the returned array).
##'
##' @return An array with 3 dimensions (state x particle x time) or 4
##'   dimensions (state x particle x group x time) for a grouped
##'   system.
##'
##' @export
dust_system_simulate <- function(sys, times, index = NULL) {
  check_is_dust_system(sys)
  ret <- sys$methods$simulate(sys$ptr, times, index, sys$grouped)
  if (!is.null(index) && !is.null(names(index))) {
    rownames(ret) <- names(index)
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
  sys$methods$reorder(sys$ptr, index)
  invisible()
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
  sys$methods$compare_data(sys$ptr, data, sys$grouped)
}


##' @export
print.dust_system <- function(x, ...) {
  cli::cli_h1("<dust_system: {x$name}>")
  if (x$grouped) {
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
  invisible(x)
}


##' @export
dim.dust_system <- function(x, ...) {
  c(x$n_state, x$n_particles, if (x$grouped) x$n_groups)
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
