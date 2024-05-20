dust_model <- function(name, env) {
  prefix <- sprintf("dust2_cpu_%s", name)
  methods <- list(name = name,
                  alloc = env[[sprintf("%s_alloc", prefix)]],
                  state = env[[sprintf("%s_state", prefix)]])
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
##' @param model
##'
##' @return An array of model state.  If your model is ungrouped, then
##'   this has two dimensions (state, particle).  If grouped, this has
##'   three dimensions (state, particle, group)
##'
##' @export
dust_model_state <- function(model) {
  model$methods$state(model$ptr, model$grouped)
}


##' @export
print.dust_model <- function(x, ...) {
  cli::cli_
  invisible(x)
}
