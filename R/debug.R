##' Control debugging of dust models.  This help page documents two
##' funtions that can be used to control if and how the debugger is
##' enabled (`dust_debug_enabled` and `dust_debug_verbosity`).  You
##' can't enter the debugger from any of these functions; it is only
##' enabled if present in your C++ code (or if using `odin2` if you
##' have enabled it).
##'
##' dust2 includes an extremely simple debugging system, and if you
##' are reading this message, there's a good chance you are inside it.
##' It is built on top of R's [browser()] and so all the usual tips,
##' tricks and issues for working with this apply.  We recommend
##' setting the R option `browserNLdisabled = TRUE` to avoid surprises
##' from presssing `<enter>`.
##'
##' * You can press `n` or `c` to proceed to the next enabled iteration
##' * You can press `Q` to quit the browser (this will end up as an
##'   error by the time you have control back)
##'
##' These commands are established by `browser`, and can't be
##' disabled.  This means that if you have an variable called `n` you
##' will need to work with it as `(n)` (i.e. in parentheses).  This
##' applies to all of browser's command variables (`c`, `f`, `n`, `s`,
##' `r` and `Q`); please see [browser] for more information.
##'
##' By the time the environment has been created, some variables from
##' your model will have been copied into the environment; you can see
##' these by running `ls()` and write expressions involving these
##' objects.  Changes that you make in R are *not* (currently)
##' propagated back into the running system.
##'
##' If you enable the debugger, you may have very many iterations to
##' get through before control is returned back to the console.  You
##' can run `dust_debug_enabled(FALSE)` to prevent entry into the
##' debugger.  We might change this interface in future.
##'
##' @title The dust debugger
##'
##' @return Both `dust_debug_enabled` and `dust_debug_verbosity`
##'   return the previous value of the option they are setting.
##'
##' @param value Logical, `TRUE` for where the debugger should be
##'   enabled, `FALSE` otherwise.
##'
##' @export
##' @rdname dust_debug
dust_debug_enabled <- function(value = TRUE) {
  if (!is.null(value)) {
    assert_scalar_logical(value)
  }
  options(dust.debug_enabled = value)
}


##' @rdname dust_debug
##' @export
##'
##' @param level The verbosity level, as a string.  This must be one
##'   of the values `quiet` (prevents informational messages),
##'   `normal` (prints a single line on entry) and `verbose` (prints
##'   several informational messages on entry).  The default is `normal`.
dust_debug_verbosity <- function(level) {
  if (!is.null(level)) {
    level <- match_value(level, c("quiet", "normal", "verbose"))
  }
  options(dust.debug_verbosity = level)
}


debug_welcome_message <- function(env) {
  verbosity <- getOption("dust.debug_verbosity", "normal")
  if (verbosity == "normal") {
    cli::cli_alert_info("dust debug: see {.help dust_debug} for help")
  } else if (verbosity == "verbose") {
    cli::cli_alert_info("dust debug")
    nms <- ls(env, all = TRUE)
    cli::cli_alert_info("{length(nms)} variable{?s} available: {nms}")
    cli::cli_alert_info("'c' to continue, 'Q' to exit")
    cli::cli_alert_info("See {.help dust_debug} for help")
  }
}


debug_env <- function(env) {
  enabled <- getOption("dust.debug_enabled", TRUE)
  ## It looks like we can use sys.frames() here to find out the
  ## environment where we started the calls from; we could use that at
  ## some point to implement something where we can prevent debugging
  ## for the rest of the calls in the current run by adding a deferred
  ## handler to that environment.  An extension for later I think.
  if (enabled) {
    debug_welcome_message(env)
    with(env, browser(skipCalls = 1L))
  }
}
