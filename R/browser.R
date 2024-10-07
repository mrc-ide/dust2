##' Control browser-based of dust models.  This help page documents
##' three funtions that can be used to control if and how the browser
##' is enabled.  You can't enter the debugger from any of these
##' functions; it is only enabled if present in your C++ code (or if
##' using `odin2` if you have enabled it).
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
##' can run `dust_debug_continue()` to prevent entry into the debugger
##' until control is passed back to the you; this means the time
##' series will run to completion and then the next time you run the
##' system the debugger will be triggered again.  Alternatively, you
##' can run `dust_debug_enabled(FALSE)` to disable all calls to the
##' debugger.
##'
##' @title The dust debugger
##'
##' @return Both `dust_browser_enabled` and `dust_browser_verbosity`
##'   return the previous value of the option they are setting.
##'
##' @param value Logical, `TRUE` for where the debugger should be
##'   enabled, `FALSE` otherwise.
##'
##' @export
##' @rdname dust_browser
##' @name dust_browser
dust_browser_enabled <- function(value = TRUE) {
  if (!is.null(value)) {
    assert_scalar_logical(value)
  }
  options(dust.browser_enabled = value)
}


##' @rdname dust_browser
##' @export
##'
##' @param level The verbosity level, as a string.  This must be one
##'   of the values `quiet` (prevents informational messages),
##'   `normal` (prints a single line on entry) and `verbose` (prints
##'   several informational messages on entry).  The default is `normal`.
dust_browser_verbosity <- function(level) {
  if (!is.null(level)) {
    level <- match_value(level, c("quiet", "normal", "verbose"))
  }
  options(dust.browser_verbosity = level)
}


##' @rdname dust_browser
##' @export
dust_browser_continue <- function() {
  parent <- browser_find_parent_env(sys.frames())
  if (is.null(parent) || environmentIsLocked(parent)) {
    cli::cli_abort(
      "Called 'dust_browser_continue()' from outside of a dust browser context")
  }
  parent$.dust_browser_continue <- TRUE
}


browser_welcome_message <- function(env, phase, time) {
  verbosity <- getOption("dust.browser_verbosity", "normal")
  if (verbosity == "normal") {
    cli::cli_alert_info(
      "dust debug ('{phase}'; time = {time}): see {.help dust_browser} for help")
  } else if (verbosity == "verbose") {
    cli::cli_alert_info("dust debug ('{phase}'; time = {time})")
    nms <- ls(env)
    cli::cli_alert_info("{length(nms)} variable{?s} available: {nms}")
    cli::cli_alert_info("'c' to continue, 'Q' to exit")
    cli::cli_alert_info(
      "Run {.run dust2::dust_browser_continue()} to stop debugging")
    cli::cli_alert_info("See {.help dust_browser} for help")
  }
}


browser_env <- function(env, phase, time) {
  enabled <- getOption("dust.browser_enabled", TRUE)
  if (enabled) {
    parent <- browser_find_parent_env(sys.frames())
    if (is.null(parent) || !isTRUE(parent$.dust_browser_continue)) {
      browser_welcome_message(env, phase, time)
      with(env, browser())
    }
  }
}


browser_find_parent_env <- function(frames, drop = 1) {
  ## Finds the outermost call to dust2, dropping the last frame(s),
  ## which we expect to be dust calls.
  frames <- drop_last(frames, drop)
  i <- vcapply(frames, function(x) environmentName(parent.env(x))) == "dust2"
  if (any(i)) {
    frames[[which(i)[[1]]]]
  } else {
    NULL
  }
}
