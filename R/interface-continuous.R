##' Create a control object for controlling the adaptive stepper for
##' systems of ordinary differential equations (ODEs). The returned
##' object can be passed into a continuous-time dust model on
##' initialisation.
##'
##' @title Create a dust_ode_control object.
##'
##' @param max_steps Maximum number of steps to take. If the
##'   integration attempts to take more steps that this, it will
##'   throw an error, stopping the integration.
##'
##' @param rtol The per-step relative tolerance.  The total accuracy
##'   will be less than this.
##'
##' @param atol The per-step absolute tolerance.
##'
##' @param step_size_min The minimum step size.  The actual minimum
##'   used will be the largest of the absolute value of this
##'   `step_size_min` or `.Machine$double.eps` (or the
##'   single-precision equivalent once we support `float`-based
##'   models).  If the integration attempts to make a step smaller than
##'   this, it will throw an error, stopping the integration.
##'
##' @param step_size_max The largest step size.  By default there is
##'   no maximum step size (`Inf`) so the solver can take as large a
##'   step as it wants to.  If you have short-lived fluctuations in
##'   your rhs that the solver may skip over by accident, then specify
##'   a smaller maximum step size here.
##'
##' @param debug_record_step_times Logical, indicating if the step
##'   times should be recorded.  This should only be enabled for
##'   debugging.  Step times can be retrieved via
##'   [dust_system_internals()].
##'
##' @param save_history Optional vector of states for which we should
##'   save index.  This is intended only for debugging.
##'
##' @export
##'
##' @return A named list of class "dust_ode_control".  Do not modify
##'   this after creation.
dust_ode_control <- function(max_steps = 10000, atol = 1e-6, rtol = 1e-6,
                             step_size_min = 0, step_size_max = Inf,
                             debug_record_step_times = FALSE,
                             save_history = NULL) {
  ctl <- list(
    max_steps = assert_scalar_size(
      max_steps, allow_zero = FALSE),
    atol = assert_scalar_positive_numeric(
      atol, allow_zero = FALSE),
    rtol = assert_scalar_positive_numeric(
      rtol, allow_zero = FALSE),
    step_size_min = assert_scalar_positive_numeric(
      step_size_min, allow_zero = TRUE),
    step_size_max = assert_scalar_positive_numeric(
      step_size_max, allow_zero = TRUE),
    save_history = assert_scalar_logical(
      save_history),
    debug_record_step_times = assert_scalar_logical(
      debug_record_step_times))
  class(ctl) <- "dust_ode_control"
  ctl
}
