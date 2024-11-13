##' Return information about OpenMP support for this machine.
##'
##' @title Information about OpenMP support
##'
##' @param check_compile Logical, indicating if we should check if we
##'   can compile an OpenMP program - this is slow the first time.
##'
##' @seealso [dust_openmp_threads()] for setting a polite number of
##'   threads.
##'
##' @return A list with information about the OpenMP support on your
##'   machine.
##'
##' * The first few elements come from the OpenMP library directly:
##'   `num_proc`, `max_threads`, `thread_limit`; these correspond to a
##'   call to the function `omp_get_<name>()` in C and
##'   `openmp_version` which is the value of the `_OPENMP` macro.
##' * A logical `has_openmp` which is `TRUE` if it looks like runtime
##'   OpenMP support is available
##' * The next elements tell you about different sources that might
##'   control the number of threads allowed to run: `mc.cores` (from
##'   the R option with the same name), `OMP_THREAD_LIMIT`,
##'   `OMP_NUM_THREADS`, `MC_CORES` (from environment variables),
##'   `limit_r` (limit computed against R-related control variables),
##'   `limit_openmp` (limit computed against OpenMP-related variables)
##'   and `limit` the smaller of `limit_r` and `limit_openmp`
##' * Finally, if you specified `check_compile = TRUE`, the logical
##'   `has_openmp_compiler` will indicate if it looks like we can
##'   compile with OpenMP.
##'
##' @export
##' @examples
##' dust_openmp_support()
dust_openmp_support <- function(check_compile = FALSE) {
  info <- openmp_info()
  if (check_compile) {
    info$has_openmp_compiler <- has_openmp_compiler()
  }
  info
}


##' Politely select a number of threads to use. See Details for the
##' algorithm used.
##'
##' There are two limits and we will take the smaller of the two.
##'
##' The first limit comes from piggy-backing off of R's normal
##' parallel configuration; we will use the `MC_CORES` environment
##' variable and `mc.cores` option as a guide to how many cores you
##' are happy to use. We take `mc.cores` first, then `MC_CORES`, which
##' is the same behaviour as `parallel::mclapply` and friends.
##'
##' The second limit comes from OpenMP. If you do not have OpenMP
##' support, then we use one thread (higher numbers have no effect at
##' all in this case). If you do have OpenMP support, we take the
##' smallest of the number of "processors" (reported by
##' `omp_get_num_procs()`) the "max threads" (reported by
##' `omp_get_max_threads()` and "thread_limit" (reported by
##' `omp_get_thread_limit()`.
##'
##' See [dust_openmp_support()] for the values of all the values that
##' go into this calculation.
##'
##' @title Select number of threads
##'
##' @param n Either `NULL` (select automatically) or an integer as
##'   your proposed number of threads.
##'
##' @param action An action to perform if `n` exceeds the maximum
##'   number of threads you can use. Options are "error" (the default,
##'   throw an error) or "fix" (print a message and reduce `n` down to
##'   the limit).
##'
##' @return An integer, indicating the number of threads that you can use
##' @export
##' @examples
##' # Default number of threads; tries to pick something friendly,
##' # erring on the conservative side.
##' dust_openmp_threads(NULL)
##'
##' # Try to pick something silly and it will be reduced for you
##' dust_openmp_threads(1000, action = "fix")
dust_openmp_threads <- function(n = NULL, action = "error") {
  info <- openmp_info()
  if (is.null(n)) {
    n <- info$limit
  } else  {
    n <- openmp_check_limit(n, info$limit, action)
  }
  n
}


has_openmp_compiler <- function() {
  if (is.null(cache$has_openmp_compiler)) {
    cache$has_openmp_compiler <- has_openmp_compiler_test()
  }
  cache$has_openmp_compiler
}


## This test uses the 'parallel' example, which as its update() method
## returns the thread number by running omp_get_thread_num()
has_openmp_compiler_test <- function() {
  workdir <- tempfile("dust_")
  dir_create(workdir)
  dir_create(file.path(workdir, "src"))
  data <- list(package = "dustopenmp",
               linking_to = "cpp11, dust2, monty",
               compiler_options = "",
               system_requirements = "R (>= 4.0.0)")
  writeLines(substitute_dust_template(data, "DESCRIPTION"),
             file.path(workdir, "DESCRIPTION"))
  writeLines(substitute_dust_template(data, "NAMESPACE"),
             file.path(workdir, "NAMESPACE"))
  writeLines(substitute_dust_template(data, "Makevars"),
             file.path(workdir, "src", "Makevars"))
  stopifnot(file.copy(dust2_file("openmp.cpp"),
                      file.path(workdir, "src"),
                      overwrite = TRUE))
  tryCatch({
    pkgbuild::compile_dll(workdir, compile_attributes = TRUE,
                          quiet = TRUE, debug = FALSE)
    env <- load_temporary_package(workdir, data$package, TRUE)
    env$openmp_get_thread_id() >= 0
  }, error = function(e) FALSE)
}


## NOTE: This does not return if the *compiler* supports OpenMP, just
## the runtime.  While we are testing that will be the same thing, but
## after installation from binary this requires really a compile time
## test of a simple OpenMP program.
openmp_info <- function() {
  env <- Sys.getenv(c("OMP_THREAD_LIMIT", "OMP_NUM_THREADS", "MC_CORES"))
  env <- set_names(as.list(as.integer(env)), names(env))
  info <- cpp_openmp_info()
  info[["mc.cores"]] <- getOption("mc.cores", NA_integer_)

  limit <- list()
  limit$limit_r <- getOption("mc.cores", as.integer(Sys.getenv("MC_CORES", 1)))
  limit$limit_openmp <- min(info$num_procs,
                            info$num_threads,
                            info$thread_limit)
  if (!info$has_openmp) {
    limit$limit_openmp <- 1L
  }
  limit$limit <- min(limit$limit_r, limit$limit_openmp)

  c(info, env, limit)
}


openmp_check_limit <- function(n, limit, action, call = parent.frame()) {
  match_value(action, c("error", "fix"))
  if (n > limit) {
    msg <- "Requested number of threads '{n}' exceeds a limit of '{limit}'"
    hint <- "See {.help [dust_openmp_threads()](dust2::dust_openmp_threads)} for details"
    if (action == "error") {
      cli::cli_abort(c(msg, i = hint), call = call)
    } else {
      cli::cli_alert_warning(msg)
      cli::cli_alert_info(hint)
      n <- limit
    }
  }
  n
}


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
