##' Create a [mcstate_model] from a [dust_filter] object.
##'
##' # Random number streams
##'
##' The short version: the seed argument that you may have passed to
##' [dust_filter_create] will be ignored when using a
##' `dust_filter_mcstate` model with algorithms from mcstate2.  You
##' should generally not worry about this, it is expected.
##'
##' When you initialise a filter, you provide a random number seed
##' argument.  Your filter will use `n_groups * (n_particles + 1)`
##' streams (one for the filter for each group, then for each group
##' one per particle).  If you run the filter directly (e.g., with
##' [dust_filter_run] then you will advance the state of the filter).
##' However, if you use the filter with mcstate2 (which is why you're
##' using `dust_filter_mcstate`) we will ignore this seeding.
##'
##' When running mcmc with `n_chains` chains, we need `n_chains *
##' n_groups * (n_particles + 1)` random number streams - that is
##' enough streams for every chain to have a filter with its own set
##' of chains.  mcstate2 will look after this for us, but the upshot
##' is that the random number state that you may have previously set
##' when building the filter will be ignored as we need to create a
##' series of suitable seeds.
##'
##' The seeds provided by mcstate will start at some point in the RNG
##' state space (2^256 possible states by default).  In an MCMC, each
##' chain will have a seed that is created by performing a
##' "long jump", moving 2^192 steps along the chain.  Then within each
##' chain we will take a series of "jumps" (2^128 steps) for each of
##' the streams we need across the groups, filters and particles.
##' This ensures indepdence across the stochastic components of the
##' system but also the reproducibility and predictability of the
##' system.  The intiial seeding performed by mcstate2 will respond to
##' R's RNG (i.e., it will follow `set.seed`) if an explicit seed is
##' not given.
##'
##' @title Create mcstate model
##'
##' @param filter A `dust_filter` object, created from [dust_filter]
##'
##' @param packer A parameter packer, which will convert between an
##'   unstructured vector of parameters as used in an MCMC into the
##'   list of parameters that the dust system requires.
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
  properties <- mcstate2::mcstate_model_properties(
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
      ## We need to expand state here as mcstate will only hand us the
      ## initial seed (i.e., the first 32 bytes of state) and we have
      ## to create the other streams by jumping.
      n_groups <- max(1, filter$n_groups)
      state <- filter_rng_state(filter$n_particles, n_groups, state)
      dust_filter_set_rng_state(filter, state)
    }
  } else {
    get_rng_state <- NULL
    set_rng_state <- NULL
  }

  ## I think that reconstructing the filter here to decouple any
  ## linkage with the provided filter would be easier to understand,
  ## so that we really copy the filter.  This just requires a copy
  ## method, which feels like it would be fairly easy to implement.
  mcstate2::mcstate_model(
    list(filter = filter,
         density = density,
         direct_sample = NULL,
         parameters = packer$parameters,
         domain = domain,
         get_rng_state = get_rng_state,
         set_rng_state = set_rng_state),
    properties)
}
