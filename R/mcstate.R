##' Create a [mcstate_model] from a [dust_filter] object.
##'
##' @title Create mcstate model
##'
##' @param filter A `dust_filter` object, created from [dust_filter]
##'
##' @param transformer A parameter transformer (unpacker?), which will
##'   convert between an unstructured vector of parameters as used in
##'   an MCMC into the list of parameters that the dust system
##'   requires.
##'
##' @return A [mcstate::mcstate_model] object
##' 
##' @export
dust_filter_mcstate <- function(filter, transformer) {
  check_is_dust_filter(filter)
  assert_is(transformer, "mcstate_transformer")

  ## We configure saving trajectories on creation I think, which then
  ## affects density.  Start without trajectories?  Realistically
  ## people don't change this much afterwards?

  ## TODO: transform needs to know about domain, I think, or at least
  ## accept information on it.

  is_stochastic <- !filter$deterministic

  density <- function(x) {
    pars <- transformer$transform(x)
    ## TODO: initial state handling here, this requires some thought,
    ## especially once we support arbitrary initial conditions, and
    ## once these come from the parameters, and once we support
    ## differentiation.  For now we require that the initialisation
    ## can come from the model itself and later we will allow
    ## overriding that - it's possible that is something that should
    ## be done at the filter level tbh.
    initial <- NULL
    dust_filter_run(filter, pars, initial)
  }

  if (is_stochastic) {
    get_rng_state <- function() {
      dust_filter_rng_state(filter)
    }
    set_rng_state <- function(state) {
      dust_filter_set_rng_state(filter, state)
    }
  } else {
    get_rng_state <- NULL
    set_rng_state <- NULL
  }

  properties <- mcstate_model_properties(
    is_stochastic = is_stochastic,
    has_direct_sample = FALSE,
    has_gradient = FALSE,
    ## We need to get this working with the transformation!
    has_parameter_groups = FALSE)
  
  mcstate2::mcstate_model(
    list(filter = filter,
         density = density,
         direct_sample = NULL,
         parameters = parameters,
         domain = domain,
         set_rng_state = NULL,
         get_rng_state = NULL),
    properties)
}
