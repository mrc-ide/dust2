## This is a basic reimplementation of a particle filter for the SIR
## model; once we have a generic model interface we can make this more
## generic.  It just redoes the same logic as in the C++ code but is
## easier to read (and quite a lot slower due to churn in state).
sir_filter_manual <- function(pars, time_start, time, dt, data, n_particles,
                              seed) {
  r <- mcstate2::mcstate_rng$new(n_streams = 1, seed = seed)
  seed <- mcstate2::mcstate_rng$new(n_streams = 1, seed = seed)$jump()$state()

  obj <- dust_model_create(sir(), pars, n_particles,
                           time = time_start, dt = dt, seed = seed)
  n_steps <- round((time - c(time_start, time[-length(time)])) / dt)
  n_state <- nrow(dust_model_state(obj))
  n_time <- length(time)

  function(pars, initial = NULL, save_history = FALSE) {
    if (!is.null(pars)) {
      dust_model_update_pars(obj, pars)
    }
    dust_model_set_time(obj, time_start)
    if (is.null(initial)) {
      dust_model_set_state_initial(obj)
    } else {
      dust_model_set_state(obj, initial)
    }
    ll <- 0
    history <- array(NA_real_, c(n_state, n_particles, n_time))
    for (i in seq_along(time)) {
      dust_model_run_steps(obj, n_steps[[i]])
      tmp <- dust_model_compare_data(obj, data[[i]])
      w <- exp(tmp - max(tmp))
      ll <- ll + log(mean(w)) + max(tmp)
      u <- r$random_real(1)
      k <- test_resample_weight(w, u) + 1L
      state <- dust_model_state(obj)
      if (save_history) {
        history[, , i] <- state
        history <- history[, k, , drop = FALSE]
      }
      dust_model_set_state(obj, state[, k])
    }
    list(log_likelihood = ll, history = if (save_history) history)
  }
}
