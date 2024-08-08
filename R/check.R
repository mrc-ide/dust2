## Distinct from the code in util_assert.R which are our generic
## validation functions, these are bits of validation that group
## together to validate inputs to dust systems.
check_pars <- function(pars, n_groups, preserve_group_dimension,
                       name = deparse(substitute(pars)), call = NULL) {
  if (preserve_group_dimension) {
    if (length(pars) != n_groups) {
      cli::cli_abort(
        paste("Expected 'pars' to have length {n_groups} to match 'n_groups',",
              "but it had length {length(pars)}"),
        arg = "pars", call = call)
    }
  } else {
    pars <- list(pars)
  }
  pars
}
