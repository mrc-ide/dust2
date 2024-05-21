dust_model <- function(name, env) {
  prefix <- sprintf("dust2_cpu_%s", name)
  ## I don't love that this requires running through sprintf() each
  ## time we create a model, but using a function for the model (see
  ## sir()), rather than an object, means that it's easier to think
  ## about the dependencies among packages.  This is also essentially
  ## how DBI works.

  methods_nms <- c("alloc",
                   "state", "set_state", "set_state_initial")

  methods <- lapply(sprintf("dust2_cpu_%s_%s", name, methods_nms),
                    function(x) env[[x]])
  names(methods) <- methods_nms
  ret <- list(name = name,
              methods = methods)
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
##' @param model A model generator object, with class
##'   `dust_model_generator`
##'
##' @param pars A list of parameters.  The format of this will depend
##'   on the model.  If `n_groups` is 1 or more, then this must be a
##'   list where each element of length `n_groups` where each element
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
dust_model_create <- function(model, pars, n_particles, n_groups = 0,
                              time = 0, dt = 1,
                              seed = NULL, deterministic = FALSE) {
  if (!inherits(model, "dust_model_generator")) {
    ## TODO: we can detect if the user has passed 'sir' rather than
    ## 'sir()' probably.
    cli::cli_abort("Expected 'model' to be a 'dust_model_generator' object")
  }
  res <- model$methods$alloc(pars, time, dt, n_particles, n_groups,
                             seed, deterministic)
  ## Here, we augment things slightly
  res$name <- model$name
  res$methods <- model$methods
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

##' Set model state from a models initial conditions.  This may depend
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


##' @export
print.dust_model_generator <- function(x, ...) {
  cli::cli_h1("<dust_model_generator: {x$name}>")
  ## Later, we might print some additional capabilities of the model
  ## here, such as if it can be used with a filter, a summary of its
  ## parameters (once we know how to access that), etc.
  cli::cli_bullets(
    i = "Use 'dust2::dust_model_create()' to create a model")
  invisible(x)
}


##' @export
print.dust_model <- function(x, ...) {
  cli::cli_h1("<dust_model: {x$name}>")
  ## Later, we might print some additional capabilities of the model
  ## here, such as if it can be used with a filter, a summary of its
  ## parameters (once we know how to access that), etc.
  invisible(x)
}


##' @export
dim.dust_model <- function(x, ...) {
  c(x$n_state, x$n_particles, if (x$grouped) x$n_groups)
}


check_is_dust_model <- function(model, call = parent.frame()) {
  if (!inherits(model, "dust_model")) {
    cli::cli_abort("Expected 'model' to be a 'dust_model' object",
                   arg = "model", call = call)
  }
}
