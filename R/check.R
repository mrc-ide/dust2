## Distinct from the code in util_assert.R which are our generic
## validation functions, these are bits of validation that group
## together to validate inputs to dust systems.
check_pars <- function(pars, n_groups, index_group, preserve_group_dimension,
                       name = deparse(substitute(pars)), call = NULL) {
  if (preserve_group_dimension) {
    expected <- if (is.null(index_group)) n_groups else length(index_group)
    if (length(pars) != expected) {
      what <- if (is.null(index_group)) "n_groups" else "index_group"
      cli::cli_abort(
        paste("Expected 'pars' to have length {expected} to match '{what}',",
              "but it had length {length(pars)}"),
        arg = "pars", call = call)
    }
  } else {
    pars <- list(pars)
  }
  pars
}
