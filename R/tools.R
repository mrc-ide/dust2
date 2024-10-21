##' Unpack state.  You will see state come out of dust2 systems in
##' several places, for example [dust_system_state], but it will
##' typically be an unstructured vector with no names; this is not very
##' useful!  However, your model knows what each element, or group of
##' elements "means".  You can unpack your state from this
##' unstructured array into a named list using this function.  The
##' same idea applies to the higher-dimensional arrays that you might
##' get if your system has multiple particles, multiple parameter
##' groups or has been run for multiple time steps.
##'
##' @title Unpack state
##'
##' @param obj A `dust_system` object (from [dust_system_create]) or
##'   `dust_likelihood` object (from [dust_filter_create] or
##'   [dust_unfilter_create]).
##'
##' @param state A state vector, matrix or array.  This might have
##'   come from [dust_system_state], [dust_likelihood_last_history], or
##'   [dust_likelihood_last_state].
##'
##' @return A named list, where each element corresponds to a logical
##'   compartment.
##'
##' @seealso [monty::monty_packer], and within that especially
##'   documentation for `$unpack()`, which powers this function.
##'
##' @rdname dust_unpack
##' @export
##' @examples
##' sir <- dust_example("sir")
##' sys <- dust_system_create(sir, list(), n_particles = 10, dt = 0.25)
##' dust_system_set_state_initial(sys)
##' t <- seq(0, 100, by = 5)
##' y <- dust_system_simulate(sys, t)
##' # The result here is a 5 x 10 x 21 matrix: 5 states by 10 particles by
##' # 21 times.
##' dim(y)
##'
##' # The 10 particles and 21 times (following t) are simple enough, but
##' # what are our 5 compartments?
##'
##' # You can use dust_unpack_state() to reshape your output as a
##' # list:
##' dust_unpack_state(sys, y)
##'
##' # Here, the list is named following the compartments (S, I, R,
##' # etc) and is a 10 x 21 matrix (i.e., the remaining dimensions
##' # from y)
##'
##' # We could apply this to the final state, which converts a 5 x 10
##' # matrix of state into a 5 element list of vectors, each with
##' # length 10:
##' s <- dust_system_state(sys)
##' dim(s)
##' dust_unpack_state(sys, s)
##'
##' # If you need more control, you can use 'dust_unpack_index' to map
##' # names to positions within the state dimension of this array
##' dust_unpack_index(sys)
dust_unpack_state <- function(obj, state) {
  get_unpacker(obj)$unpack(state)
}


##' @rdname dust_unpack
##' @export
dust_unpack_index <- function(obj) {
  get_unpacker(obj)$index()
}


get_unpacker <- function(obj, call = parent.frame()) {
  if (inherits(obj, "dust_likelihood")) {
    if (is.null(obj$packer_state)) {
      cli::cli_abort(
        c("Packer is not yet ready",
          i = "Likelihood has not yet been run"))
    }
  } else if (!inherits(obj, "dust_system")) {
    cli::cli_abort(
      "Expected 'obj' to be a 'dust_system' or a 'dust_likelihood'",
      arg = "obj", call = call)
  }
  obj$packer_state
}
