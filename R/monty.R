##' Create a [monty_model] from a `dust_likelihood` object.
##'
##' # Random number streams
##'
##' This section is only relevant if your likelihood object is a
##' particle filter, and therefore uses random numbers.
##'
##' The short version: the seed argument that you may have passed to
##' [dust_filter_create] will be ignored when using a
##' `dust_likelihood_monty` model with algorithms from monty.  You
##' should generally not worry about this, it is expected.
##'
##' When you initialise a filter, you provide a random number seed
##' argument.  Your filter will use `n_groups * (n_particles + 1)`
##' streams (one for the filter for each group, then for each group
##' one per particle).  If you run the filter directly (e.g., with
##' [dust_likelihood_run] then you will advance the state of the filter).
##' However, if you use the filter with monty (which is why you're
##' using `dust_likelihood_monty`) we will ignore this seeding.
##'
##' When running MCMC with `n_chains` chains, we need `n_chains *
##' n_groups * (n_particles + 1)` random number streams - that is
##' enough streams for every chain to have a filter with its own set
##' of chains.  monty will look after this for us, but the upshot
##' is that the random number state that you may have previously set
##' when building the filter will be ignored as we need to create a
##' series of suitable seeds.
##'
##' The seeds provided by monty will start at some point in the RNG
##' state space (2^256 possible states by default).  In an MCMC, each
##' chain will have a seed that is created by performing a
##' "long jump", moving 2^192 steps along the chain.  Then within each
##' chain we will take a series of "jumps" (2^128 steps) for each of
##' the streams we need across the groups, filters and particles.
##' This ensures independence across the stochastic components of the
##' system but also the reproducibility and predictability of the
##' system.  The initial seeding performed by monty will respond to
##' R's RNG (i.e., it will follow `set.seed`) if an explicit seed is
##' not given.
##'
##' @title Create monty model
##'
##' @param obj A `dust_likelihood` object, created from
##'   [dust_filter_create] or [dust_unfilter_create]
##'
##' @param packer A parameter packer, which will convert between an
##'   unstructured vector of parameters as used in an MCMC into the
##'   list of parameters that the dust system requires.
##'
##' @param initial Optionally, a function to create initial conditions
##'   from unpacked parameters.
##'
##' @param domain Optionally, domain information to pass into the
##'   model.  If given, this is a two column matrix with row names
##'   corresponding to the parameter names in `packer`, the first
##'   column representing the lower bound and the second column
##'   representing the upper bound.  You do not need to specify
##'   parameters that have a domain of `(-Inf, Inf)` as this is
##'   assumed.  We use [monty::monty_domain_expand] to expand
##'   logical parameters, so if you have a vector-valued parameter `b`
##'   and a domain entry called `b` we will expand this to all
##'   elements of `b`.
##'
##' @param failure_is_impossible Logical, indicating if an error while
##'   computing the likelihood should be treated as a log-density of
##'   `-Inf` (i.e., that this point is impossible).  This is a big
##'   hammer to use, and you would be better off using the domain
##'   (with reflecting boundaries) or the priors to control this if
##'   possible.  However, sometimes you can have integration failures
##'   with very high parameter values, or just other pathological
##'   parameter sets where, once you understand the model, giving up
##'   on that parameter set and continuing is the best option.
##'
##' @param save_state Logical, indicating if the state should be saved
##'   at the time series.
##'
##' @param save_trajectories Logical, indicating if particle trajectories
##'   should be saved.  This can also be a character vector indicating
##'   the logical compartments, or an integer vector being an index
##'   into the state.
##'
##' @return A [monty::monty_model] object
##'
##' @export
dust_likelihood_monty <- function(obj, packer, initial = NULL, domain = NULL,
                                  failure_is_impossible = FALSE,
                                  save_state = FALSE,
                                  save_trajectories = FALSE) {
  assert_is(obj, "dust_likelihood")

  is_grouped <- !is.null(obj$groups)
  if (is_grouped) {
    assert_is(packer, "monty_packer_grouped")
    if (!identical(packer$groups(), obj$groups)) {
      cli::cli_abort(
        c("Groups for 'packer' do not match those of 'obj'",
          i = "'obj' has: {squote(obj$groups)}",
          x = "'packer' has: {squote(packer$groups())}"),
        arg = "packer")
    }
  } else {
    assert_is(packer, "monty_packer")
  }

  domain <- monty::monty_domain_expand(domain, packer)
  save_trajectories <- validate_save_trajectories(save_trajectories)
  observer <- dust_observer(obj, save_state, save_trajectories$enabled)

  ## We configure saving trajectories on creation I think, which then
  ## affects density.  Start without trajectories?  Realistically
  ## people don't change this much afterwards?  This will be easiest
  ## to think about once we sort out how models know about their variables.

  ## Supporting parameter groups requires some effort with packers,
  ## movement there will come from monty, and probably we'll end up
  ## subclassing the packer.
  properties <- monty::monty_model_properties(
    is_stochastic = !obj$deterministic,
    has_direct_sample = FALSE,
    has_gradient = obj$deterministic && obj$has_adjoint,
    allow_multiple_parameters = FALSE,
    has_observer = !is.null(observer),
    has_parameter_groups = FALSE)

  gradient <- NULL

  env <- environment()

  if (obj$deterministic) {
    ## We might move this cache into the unfilter itself, but that
    ## will take a bit of an effort to get right.  The other option
    ## will be a `dust_likelihood_gradient` function but that then
    ## requires that we have another one for getting the last density
    ## back out!
    ##
    ## The main issue here is if two different bits of code are
    ## manipulating the object, then we might end up with different
    ## parameters in the cache, so rather than using an environment
    ## for the cache, we store them against the object's ptr field,
    ## which ensures that the internal structure is always in sync -
    ## we can change this later with no obvious external effect, noone
    ## but us should be depending on the internal structures of
    ## object.
    density <- function(x) {
      pars <- packer_unpack(packer, x)
      was_uninitialised <- dust_likelihood_ensure_initialised(obj, pars)
      needs_trajectories_index <- was_uninitialised &&
        save_trajectories$enabled && !is.null(save_trajectories$subset)
      if (needs_trajectories_index) {
        env$save_trajectories$index <-
          obj$packer_state$subset(save_trajectories$subset)$index
      }
      if (!identical(x, attr(obj$ptr, "last_pars"))) {
        ret <- dust_likelihood_run(
          obj,
          pars,
          initial = if (is.null(initial)) NULL else initial(pars),
          save_trajectories = save_trajectories$enabled,
          index_state = save_trajectories$index)
        attr(obj$ptr, "last_density") <- if (is_grouped) sum(ret) else ret
        attr(obj$ptr, "last_gradient") <- NULL
      }
      attr(obj$ptr, "last_density")
    }

    if (properties$has_gradient) {
      gradient <- function(x) {
        density(x)
        if (is.null(attr(obj$ptr, "last_gradient"))) {
          attr(obj$ptr, "last_gradient") <-
            dust_likelihood_last_gradient(obj)
        }
        attr(obj$ptr, "last_gradient")
      }
    }
  } else {
    density <- function(x) {
      pars <- packer_unpack(packer, x)
      was_uninitialised <- dust_likelihood_ensure_initialised(obj, pars)
      needs_trajectories_index <- was_uninitialised &&
        save_trajectories$enabled && !is.null(save_trajectories$subset)
      if (needs_trajectories_index) {
        env$save_trajectories$index <-
          obj$packer_state$subset(save_trajectories$subset)$index
      }
      ll <- dust_likelihood_run(
        obj,
        pars,
        initial = if (is.null(initial)) NULL else initial(pars),
        save_trajectories = save_trajectories$enabled,
        index_state = save_trajectories$index)
      if (is_grouped) sum(ll) else ll
    }
  }

  if (failure_is_impossible) {
    density <- protect(density, -Inf)
  }

  if (properties$is_stochastic) {
    get_rng_state <- function() {
      dust_likelihood_rng_state(obj)
    }
    set_rng_state <- function(state) {
      ## We need to expand state here as monty will only hand us the
      ## initial seed (i.e., the first 32 bytes of state) and we have
      ## to create the other streams by jumping.
      n_groups <- max(1, obj$n_groups)
      state <- filter_rng_state(obj$n_particles, n_groups, state)
      dust_likelihood_set_rng_state(obj, state)
    }
  } else {
    get_rng_state <- NULL
    set_rng_state <- NULL
  }

  ## I think that reconstructing the filter here to decouple any
  ## linkage with the provided filter would be easier to understand,
  ## so that we really copy the filter.  This just requires a copy
  ## method, which feels like it would be fairly easy to implement.
  ##
  ## TODO: we have a copy method now, so is this comment ready to deal
  ## with?
  monty::monty_model(
    list(likelihood = obj,
         density = density,
         direct_sample = NULL,
         gradient = gradient,
         parameters = packer$names(),
         domain = domain,
         observer = observer,
         get_rng_state = get_rng_state,
         set_rng_state = set_rng_state),
    properties)
}


packer_unpack <- function(packer, x) {
  if (is.matrix(x)) {
    lapply(seq_len(ncol(x)), function(i) packer$unpack(x[, i]))
  } else {
    packer$unpack(x)
  }
}


validate_save_trajectories <- function(save_trajectories,
                                       call = parent.frame()) {
  if (isFALSE(save_trajectories)) {
    return(list(enabled = FALSE))
  }

  if (isTRUE(save_trajectories)) {
    return(list(enabled = TRUE, index = NULL, subset = NULL))
  }

  if (!is.character(save_trajectories)) {
    cli::cli_abort(
      paste("Invalid value for 'save_trajectories', expected a scalar",
            "logical or a character vector"),
      arg = "save_trajectories", call = call)
  }

  list(enabled = TRUE, index = NULL, subset = save_trajectories)
}


dust_observer <- function(obj, save_state, save_trajectories,
                          call = parent.frame()) {
  assert_scalar_logical(save_state, call = call)

  if (!save_trajectories && !save_state) {
    return(NULL)
  }

  observe <- function() {
    ret <- list()
    if (save_state) {
      ret$state <- dust_likelihood_last_state(
        obj, select_random_particle = TRUE)
    }

    if (save_trajectories) {
      ret$trajectories <- dust_likelihood_last_trajectories(
        obj,
        select_random_particle = TRUE)
    }

    ret
  }

  monty::monty_observer(observe)
}
