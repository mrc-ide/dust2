## This is a basic reimplementation of a particle filter for the SIR
## system; once we have a generic system interface we can make this more
## generic.  It just redoes the same logic as in the C++ code but is
## easier to read (and quite a lot slower due to churn in state).
sir_filter_manual <- function(pars, time_start, data, dt, n_particles,
                              seed) {
  r <- mcstate2::mcstate_rng$new(n_streams = 1, seed = seed)
  seed <- mcstate2::mcstate_rng$new(n_streams = 1, seed = seed)$jump()$state()

  obj <- dust_system_create(sir(), pars, n_particles,
                            time = time_start, dt = dt, seed = seed)
  n_state <- nrow(dust_system_state(obj))
  n_time <- nrow(data)
  time <- data$time

  function(pars, initial = NULL, save_history = FALSE) {
    if (!is.null(pars)) {
      dust_system_update_pars(obj, pars)
    }
    dust_system_set_time(obj, time_start)
    if (is.null(initial)) {
      dust_system_set_state_initial(obj)
    } else {
      dust_system_set_state(obj, initial)
    }
    ll <- 0
    history <- array(NA_real_, c(n_state, n_particles, n_time))
    for (i in seq_along(time)) {
      dust_system_run_to_time(obj, time[[i]])
      tmp <- dust_system_compare_data(obj, as.list(data[i, ]))
      w <- exp(tmp - max(tmp))
      ll <- ll + log(mean(w)) + max(tmp)
      u <- r$random_real(1)
      k <- test_resample_weight(w, u) + 1L
      state <- dust_system_state(obj)
      if (save_history) {
        history[, , i] <- state
        history <- history[, k, , drop = FALSE]
      }
      dust_system_set_state(obj, state[, k])
    }
    list(log_likelihood = ll, history = if (save_history) history)
  }
}


skip_for_compilation <- function() {
  skip_on_cran()
  tryCatch(
    stop_unless_installed(dust_compile_needs()),
    error = function(e) testthat::skip(conditionMessage(e)))
}


logistic_analytic <- function(r, k, times, y0) {
  len <- max(length(r), length(k), length(y0))
  vapply(times, function(t) k / (1 + (k / y0 - 1) * exp(-r * t)),
         numeric(len))
}
