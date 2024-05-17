## This is a basic reimplementation of a particle filter for the SIR
## model; once we have a generic model interface we can make this more
## generic.  It just redoes the same logic as in the C++ code but is
## easier to read (and quite a lot slower due to churn in state).
sir_filter_manual <- function(pars, time_start, time, dt, data, n_particles,
                              seed) {
  r <- mcstate2::mcstate_rng$new(n_streams = 1, seed = seed)
  seed <- mcstate2::mcstate_rng$new(n_streams = 1, seed = seed)$jump()$state()

  obj <- dust2_cpu_sir_alloc(pars, time_start, dt, n_particles, 0, seed, FALSE)
  ptr <- obj[[1]]
  n_steps <- round((time - c(time_start, time[-length(time)])) / dt)

  function(pars) {
    if (!is.null(pars)) {
      dust2_cpu_sir_update_pars(ptr, pars, FALSE)
    }
    expect_null(dust2_cpu_walk_set_time(ptr, time_start))
    dust2_cpu_sir_set_state_initial(ptr)
    ll <- 0
    for (i in seq_along(time)) {
      dust2_cpu_sir_run_steps(ptr, n_steps[[i]])
      tmp <- dust2_cpu_sir_compare_data(ptr, data[[i]], FALSE)
      w <- exp(tmp - max(tmp))
      ll <- ll + log(mean(w)) + max(tmp)
      u <- r$random_real(1)
      k <- test_resample_weight(w, u) + 1L
      state <- dust2_cpu_sir_state(ptr, FALSE)
      dust2_cpu_sir_set_state(ptr, state[, k], FALSE)
    }
    ll
  }
}
