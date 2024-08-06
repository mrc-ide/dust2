## We should have the system report back if it supports openmp at all,
## and pass that in here too, because that deserves a warning.
##
## We should also pass in the number of cores available and warn if
## the user goes over this.
##
## We should make all the warnings here tuneable, because it could get
## annoying.  They certainly should not be errors, because the user
## might have a good reason for doing it.
##
## This should do for now as a placeholder though?
check_n_threads <- function(n_threads, n_particles, n_groups, call = NULL) {
  assert_scalar_size(n_threads, allow_zero = FALSE, arg = "n_threads",
                     call = call)
  n_particles_total <- max(1, max(1, n_groups) * n_particles)
  if (n_threads > n_particles_total) {
    cli::cli_warn(
      paste("Reducing 'n_threads' from requested {n_threads} to",
            "{n_particles_total}, to match the total number of particles"))
    n_threads <- n_particles_total
  }
  as.integer(n_threads)
}
