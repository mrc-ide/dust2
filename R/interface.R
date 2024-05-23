dust_model <- function(name, env = parent.env(parent.frame())) {
  prefix <- sprintf("dust2_cpu_%s", name)
  ## I don't love that this requires running through sprintf() each
  ## time we create a model, but using a function for the model (see
  ## sir()), rather than an object, means that it's easier to think
  ## about the dependencies among packages.  This is also essentially
  ## how DBI works.
  get_methods <- function(nms, prefix) {
    set_names(
      lapply(sprintf("dust2_cpu_%s_%s", prefix, nms), function(x) env[[x]]),
      nms)
  }

  methods_core <- c("alloc",
                    "state", "set_state", "set_state_initial",
                    "time", "set_time",
                    "rng_state",
                    "update_pars",
                    "run_steps", "run_to_time", "simulate",
                    "reorder")
  methods_compare <- "compare_data"
  methods <- get_methods(c(methods_core, methods_compare), name)
  ok <- !vapply(methods[methods_core], is.null, TRUE)
  stopifnot(all(ok))

  properties <- list(
    has_compare = !is.null(methods$compare_data))

  if (properties$has_compare) {
    methods_unfilter <- c("alloc", "run", "last_history")
    methods$unfilter <-
      get_methods(methods_unfilter, sprintf("%s_unfilter", name))
    methods_filter <- c("alloc", "run", "rng_state")
    methods$filter <-
      get_methods(methods_filter, sprintf("%s_filter", name))
  }

  ret <- list(name = name,
              methods = methods,
              properties = properties)
  ## TODO: check that alloc exists, then go through and add
  ## properties.
  class(ret) <- "dust_model_generator"
  ret
}


##' Create a dust object from a model generator.  This allocates a
##' model and sets an initial set of parameters.  Once created you can
##' use other dust functions to interact with it.
##'
##' @title Create a dust object
##'
##' @param generator A model generator object, with class
##'   `dust_model_generator`
##'
##' @param pars A list of parameters.  The format of this will depend
##'   on the model.  If `n_groups` is 1 or more, then this must be a
##'   list of length `n_groups` where each element
##'   is a list of parameters for your model.
##'
##' @param n_particles The number of particles to create.
##'
##' @param n_groups Optionally, the number of parameter groups
##'
##' @param time The initial time, defaults to 0
##'
##' @param dt The time step for the model, defaults to 1
##'
##' @param seed Optionally, a seed.  Otherwise we respond to R's RNG
##'   seed on initialisation.
##'
##' @param deterministic Logical, indicating if the model should be
##'   allocated in deterministic mode.
##'
##' @return A `dust_model` object, with opaque format.
##'
##' @export
dust_model_create <- function(generator, pars, n_particles, n_groups = 0,
                              time = 0, dt = 1,
                              seed = NULL, deterministic = FALSE) {
  check_is_dust_model_generator(generator, substitute(generator))
  res <- generator$methods$alloc(pars, time, dt, n_particles, n_groups,
                                 seed, deterministic)
  ## Here, we augment things slightly
  res$name <- generator$name
  res$n_particles <- as.integer(n_particles)
  res$n_groups <- as.integer(max(n_groups), 1)
  res$deterministic <- deterministic
  res$methods <- generator$methods
  res$properties <- generator$properties
  class(res) <- "dust_model"
  res
}


##' Extract model state
##'
##' @title Extract model state
##'
##' @param model A `dust_model` object
##'
##' @return An array of model state.  If your model is ungrouped, then
##'   this has two dimensions (state, particle).  If grouped, this has
##'   three dimensions (state, particle, group)
##'
##' @seealso [dust_model_set_state()] for setting state and
##'   [dust_model_set_state_initial()] for setting state to the
##'   model-specific initial conditions.
##'
##' @export
dust_model_state <- function(model) {
  check_is_dust_model(model)
  model$methods$state(model$ptr, model$grouped)
}


##' Set model state.  Takes a multidimensional array (2- or 3d
##' depending on if the model is grouped or not).  Dimensions of
##' length 1 will be recycled as appropriate.
##'
##' @title Set model state
##'
##' @inheritParams dust_model_state
##'
##' @param state A matrix or array of state.  If ungrouped, the
##'   dimension order expected is state x particle.  If grouped the
##'   order is state x particle x group.
##'
##' @return Nothing, called for side effects only
##' @export
dust_model_set_state <- function(model, state) {
  check_is_dust_model(model)
  model$methods$set_state(model$ptr, state, model$grouped)
  invisible()
}


##' Set model state from a model's initial conditions.  This may depend
##' on the current time.
##'
##' @title Set model state to initial conditions
##'
##' @inheritParams dust_model_state
##'
##' @return Nothing, called for side effects only
##' @export
dust_model_set_state_initial <- function(model) {
  check_is_dust_model(model)
  model$methods$set_state_initial(model$ptr)
  invisible()
}


##' Fetch the current time from the model.
##'
##' @title Fetch model time
##'
##' @inheritParams dust_model_state
##'
##' @return A single numeric value
##' @seealso [dust_model_set_time]
##' @export
dust_model_time <- function(model) {
  check_is_dust_model(model)
  model$methods$time(model$ptr)
}


##' Set time into the model.  This updates the time to the provided
##' value but does not affect the state.  You may want to call
##' [dust_model_set_state] or [dust_model_set_state_initial] after
##' calling this.
##'
##' @title Set model time
##'
##' @inheritParams dust_model_time
##'
##' @param time The time to set.  Currently this must be an
##'   integer-like value, but in future we will allow setting to any
##'   multiple of `dt`.
##'
##' @return Nothing, called for side effects only
##' @export
dust_model_set_time <- function(model, time) {
  check_is_dust_model(model)
  model$methods$set_time(model$ptr, time)
  invisible()
}


##' Fetch the random number generator (RNG) state from the model.
##'
##' @title Fetch rng state
##'
##' @inheritParams dust_model_state
##'
##' @return A raw vector, this could be quite long.
##'
##' @seealso You can pass the state you get back from this function as
##'   the seed object to `dust_model_create`.  In future,
##'   `dust_model_set_rng_state` will also let you set state into an
##'   existing model.
##'
##' @export
dust_model_rng_state <- function(model) {
  check_is_dust_model(model)
  model$methods$rng_state(model$ptr)
}


##' Update parameters used by the model.  This can be used to update a
##' subset of parameters that do not change the extent of the model,
##' and will be potentially faster than creating a new model object.
##'
##' @title Update parameters
##'
##' @inheritParams dust_model_state
##'
##' @param pars Parameters to set into the model.
##'
##' @return Nothing, called for side effects only
##' @export
dust_model_update_pars <- function(model, pars) {
  check_is_dust_model(model)
  model$methods$update_pars(model$ptr, pars, model$grouped)
  invisible()
}


##' Run a model, advancing time and the state by repeatedly running
##' its `update` method.  You can advance a model either a fixed
##' (positive) number of steps, or up to a time (which must be in the
##' future).
##'
##' @title Run model
##'
##' @inheritParams dust_model_state
##'
##' @param steps The number of steps to run forward
##'
##' @return Nothing, called for side effects only
##'
##' @export
##'
##' @rdname dust_model_run
dust_model_run_steps <- function(model, steps) {
  check_is_dust_model(model)
  model$methods$run_steps(model$ptr, steps)
  invisible()
}


##' @param time Time to run to
##' @export
##' @rdname dust_model_run
dust_model_run_to_time <- function(model, time) {
  check_is_dust_model(model)
  model$methods$run_to_time(model$ptr, time)
  invisible()
}


##' Simulate a model over a series of times, returning an array of
##' output.  This output can be quite large, so you may filter states
##' according to some index.
##'
##' @title Simulate model
##'
##' @inheritParams dust_model_state
##'
##' @param times A vector of times.  They must be increasing, and the
##'   first time must be no less than the current model time
##'   (as reported by [dust_model_time])
##'
##' @param index An optional index of states to extract.  If given,
##'   then we subset the model state on return.  You can use this to
##'   return fewer model states than the model ran with, to reorder
##'   states, or to name them on exit (names present on the index will
##'   be copied into the rownames of the returned array).
##'
##' @return An array with 3 dimensions (state x particle x time) or 4
##'   dimensions (state x particle x group x time) for a grouped
##'   model.
##'
##' @export
dust_model_simulate <- function(model, times, index = NULL) {
  check_is_dust_model(model)
  ret <- model$methods$simulate(model$ptr, times, index, model$grouped)
  if (!is.null(index) && !is.null(names(index))) {
    rownames(ret) <- names(index)
  }
  ret
}


##' Reorder states within a model.  This function is primarily used
##' for debugging and may be removed from the interface if it is not
##' generally useful.
##'
##' @title Reorder states
##'
##' @inheritParams dust_model_state
##'
##' @param index The parameter ordering.  For an ungrouped model this
##'   is a vector where each element is the parameter index (if
##'   element `i` is `j` then after reordering the `i`th particle will
##'   have the state previously used by `j`).  All elements must lie
##'   in `[1, n_particles]`, repetition of an index is allowed (so
##'   that many new particles may have the state as one old particle).
##'   If the model is grouped, `index` must be a matrix with
##'   `n_particles` rows and `n_groups` columns, with each column
##'   corresponding to the reordering for a group.
##'
##' @return Nothing, called for side effects only.
##'
##' @export
dust_model_reorder <- function(model, index) {
  check_is_dust_model(model)
  model$methods$reorder(model$ptr, index)
  invisible()
}


##' Compare current model state against data.  This is only supported
##' for models that have 'compare_data' support (i.e., the model
##' definition includes a `compare_data` method).  The current state
##' in the model ([dust_model_state]) is compared against the data
##' provided as `data`.
##'
##' @title Compare model state against data
##'
##' @inheritParams dust_model_state
##'
##' @param data The data to compare against. If the model is ungrouped
##'   then `data` is a list with elements corresponding to whatever
##'   your model requires.  If your model is grouped, this should be a
##'   list with as many elements as your model has groups, with each
##'   element corresponding to the data your model requires.
##'
##' @return A numeric vector with as many elements as your model has
##'   groups, corresponding to the log-likelihood of the data for each
##'   group.
##'
##' @export
dust_model_compare_data <- function(model, data) {
  check_is_dust_model(model)
  if (!model$properties$has_compare) {
    ## This moves into something general soon?
    cli::cli_abort(
      paste("Can't compare against data; the '{model$name}' model does not",
            "have 'compare_data' support"),
      arg = "model")
  }
  model$methods$compare_data(model$ptr, data, model$grouped)
}


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
##' @inheritParams dust_model_simulate
##'
##' @return A `dust_unfilter` object, which can be used with
##'   [dust_unfilter_run]
##'
##' @export
dust_unfilter_create <- function(generator, pars, time_start, time, data,
                                 n_particles = 1, n_groups = 0,
                                 dt = 1, index = NULL) {
  check_is_dust_model_generator(generator)
  if (!generator$properties$has_compare) {
    ## This moves into something general soon?
    cli::cli_abort(
      paste("Can't create unfilter; the '{generator$name}' model does",
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
##' @param save_history Logical, indicating if the simulation history
##'   should be saved while the simulation runs; this has a small
##'   overhead in runtime and in memory.  History (particle
##'   trajectories) will be saved at each time in the filter.  If the
##'   filter was constructed using a non-`NULL` `index` parameter, the
##'   the history is restricted to these states.
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


##' Create a particle filter object
##'
##' @title Create a particle filter
##'
##' @param generator A model generator object, with class
##'   `dust_model_generator`.  The model must support `compare_data`
##'   to be used with this function.
##'
##' @param pars Initial parameters for the model.  This should be the
##'   *full* set of parameters; running the filter may update subsets
##'   using the model's underlying `update_pars` method.
##'
##' @param time_start The start time for the simulation - this is
##'   typically before the first data point.  Must be an integer-like
##'   value.
##'
##' @param time A vector of times, each of which has a corresponding
##'   entry in `data`.  The model will stop at each of these times to
##'   compute the likelihood using the compare function.
##'
##' @param data The data to compare against.  This must be a list with
##'   the same length as `time`, each element of which corresponds to
##'   the data required for the model.  If the model is ungrouped then
##'   each element of `data` is a list with elements corresponding to
##'   whatever your model requires.  If your model is grouped, this
##'   should be a list with as many elements as your model has groups,
##'   with each element corresponding to the data your model requires.
##'   We will likely introduce a friendlier data.frame based input
##'   soon.
##'
##' @param n_particles The number of particles to run.  Larger numbers
##'   will give lower variance in the likelihood estimate but run more
##'   slowly.
##'
##' @inheritParams dust_model_create
##'
##' @return A `dust_unfilter` object, which can be used with
##'   [dust_unfilter_run]
##'
##' @export
dust_filter_create <- function(generator, pars, time_start, time, data,
                               n_particles, n_groups = 0, dt = 1,
                               seed = NULL) {
  check_is_dust_model_generator(generator)
  if (!generator$properties$has_compare) {
    ## This moves into something general soon?
    cli::cli_abort(
      paste("Can't create filter; the '{generator$name}' model does",
            "not have 'compare_data' support"),
      arg = "generator")
  }
  res <- generator$methods$filter$alloc(pars, time_start, time, dt, data,
                                        n_particles, n_groups, seed)
  res$name <- generator$name
  res$n_particles <- as.integer(n_particles)
  res$n_groups <- as.integer(max(n_groups), 1)
  res$deterministic <- FALSE
  res$methods <- generator$methods$filter
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
##'   provided, the model initial conditions are used.
##'
##' @return A vector of likelihood values, with as many elements as
##'   there are groups.
##'
##' @export
dust_filter_run <- function(filter, pars = NULL, initial = NULL) {
  check_is_dust_filter(filter)
  filter$methods$run(filter$ptr, pars, initial, filter$grouped)
}


##' Get random number generator (RNG) state from the particle filter.
##'
##' @title Get filter RNG state
##'
##' @inheritParams dust_filter_run
##'
##' @return A raw vector, this could be quite long.  Later we will
##'   describe how you might reseed a filter or model with this state.
##'
##' @export
dust_filter_rng_state <- function(filter) {
  check_is_dust_filter(filter)
  filter$methods$rng_state(filter$ptr)
}


##' @export
print.dust_model_generator <- function(x, ...) {
  cli::cli_h1("<dust_model_generator: {x$name}>")
  ## Later, we might print some additional capabilities of the model
  ## here, such as if it can be used with a filter, a summary of its
  ## parameters (once we know how to access that), etc.
  cli::cli_alert_info(
    "Use 'dust2::dust_model_create()' to create a model with this generator")
  invisible(x)
}


##' @export
print.dust_model <- function(x, ...) {
  cli::cli_h1("<dust_model: {x$name}>")
  if (x$grouped) {
    cli::cli_alert_info(paste(
      "{x$n_state} state x {x$n_particles} particle{?s} x",
      "{x$n_groups} group{?s}"))
  } else {
    cli::cli_alert_info("{x$n_state} state x {x$n_particles} particle{?s}")
  }
  if (x$deterministic) {
    cli::cli_bullets(c(
      i = "This model is deterministic"))
  }
  if (x$properties$has_compare) {
    cli::cli_bullets(c(
      i = "This model has 'compare_data' support"))
  }
  ## Later, we might print some additional capabilities of the model
  ## here, such as if it can be used with a filter, a summary of its
  ## parameters (once we know how to access that), etc.
  invisible(x)
}


##' @export
dim.dust_model <- function(x, ...) {
  c(x$n_state, x$n_particles, if (x$grouped) x$n_groups)
}


check_is_dust_model_generator <- function(generator, called_as,
                                          call = parent.frame()) {
  if (!inherits(generator, "dust_model_generator")) {
    hint <- NULL
    if (is_uncalled_generator(generator) && is.symbol(called_as)) {
      hint <- c(
        i = "Did you mean '{deparse(called_as)}()' (i.e., with parentheses)")
    }
    cli::cli_abort(
      c("Expected 'generator' to be a 'dust_model_generator' object",
        hint),
      arg = "generator")
  }
}


check_is_dust_model <- function(model, call = parent.frame()) {
  if (!inherits(model, "dust_model")) {
    cli::cli_abort("Expected 'model' to be a 'dust_model' object",
                   arg = "model", call = call)
  }
}


check_is_dust_unfilter <- function(unfilter, call = parent.frame()) {
  if (!inherits(unfilter, "dust_unfilter")) {
    cli::cli_abort("Expected 'unfilter' to be a 'dust_unfilter' object",
                   arg = "unfilter", call = call)
  }
}


check_is_dust_filter <- function(filter, call = parent.frame()) {
  if (!inherits(filter, "dust_filter")) {
    cli::cli_abort("Expected 'filter' to be a 'dust_filter' object",
                   arg = "filter", call = call)
  }
}

is_uncalled_generator <- function(model) {
  if (!is.function(model)) {
    return(FALSE)
  }
  code <- body(model)
  rlang::is_call(code, "{") &&
    length(code) == 2 &&
    rlang::is_call(code[[2]], "dust_model")
}
