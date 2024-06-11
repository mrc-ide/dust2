##' Create a [mcstate_model] from a [dust_filter] object.
##'
##' @title Create mcstate model
##'
##' @param filter A `dust_filter` object, created from [dust_filter]
##'
##' @param packer A parameter packer (unpacker?), which will
##'   convert between an unstructured vector of parameters as used in
##'   an MCMC into the list of parameters that the dust system
##'   requires.
##'
##' @param initial Optionally, a function to create initial conditions
##'   from unpacked parameters.
##'
##' @return A [mcstate::mcstate_model] object
##'
##' @export
dust_filter_mcstate <- function(filter, packer, initial = NULL) {
  check_is_dust_filter(filter)
  assert_is(packer, "mcstate_packer")

  ## We configure saving trajectories on creation I think, which then
  ## affects density.  Start without trajectories?  Realistically
  ## people don't change this much afterwards?  This will be easiest
  ## to think about once we sort out how models know about their variables.

  ## TODO: Also accept information about domain to add here?  This
  ## should be in terms of the packed vector.
  domain <- NULL

  ## Supporting parameter groups requires some effort with packers,
  ## movement there will come from mcstate, and probably we'll end up
  ## subclassing the packer.
  properties <- mcstate_model_properties(
    is_stochastic = !filter$deterministic,
    has_direct_sample = FALSE,
    has_gradient = FALSE,
    has_parameter_groups = FALSE)

  density <- function(x) {
    pars <- packer$unpack(x)
    dust_filter_run(filter,
                    pars,
                    initial = if (is.null(initial)) NULL else initial(pars))
  }

  if (properties$is_stochastic) {
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

  mcstate2::mcstate_model(
    list(filter = filter,
         density = density,
         direct_sample = NULL,
         parameters = parameters,
         domain = domain,
         set_rng_state = get_rng_state,
         get_rng_state = set_rng_state),
    properties)
}
