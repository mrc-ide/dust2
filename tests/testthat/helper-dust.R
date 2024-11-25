## This is a basic reimplementation of a particle filter for the SIR
## system; once we have a generic system interface we can make this more
## generic.  It just redoes the same logic as in the C++ code but is
## easier to read (and quite a lot slower due to churn in state).
filter_manual <- function(generator, pars, time_start, data, dt, n_particles,
                          seed) {
  r <- monty::monty_rng$new(n_streams = 1, seed = seed)
  seed <- monty::monty_rng$new(n_streams = 1, seed = seed)$jump()$jump()$state()

  obj <- dust_system_create(generator, pars, n_particles,
                            time = time_start, dt = dt, seed = seed)
  n_state <- nrow(dust_system_state(obj))
  n_time <- nrow(data)
  time <- data$time

  function(pars, initial = NULL, save_trajectories = FALSE) {
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
    trajectories <- array(NA_real_, c(n_state, n_particles, n_time))
    for (i in seq_along(time)) {
      dust_system_run_to_time(obj, time[[i]])
      tmp <- dust_system_compare_data(obj, as.list(data[i, ]))
      w <- exp(tmp - max(tmp))
      ll <- ll + log(mean(w)) + max(tmp)
      u <- r$random_real(1)
      k <- test_resample_weight(w, u) + 1L
      state <- dust_system_state(obj)
      if (save_trajectories) {
        trajectories[, , i] <- state
        trajectories <- trajectories[, k, , drop = FALSE]
      }
      dust_system_set_state(obj, state[, k])
    }
    list(log_likelihood = ll,
         trajectories = if (save_trajectories) trajectories)
  }
}

sir_filter_manual <- function(...) {
  filter_manual(sir, ...)
}


malaria_filter_manual <- function(...) {
  filter_manual(malaria, ...)
}


skip_for_compilation <- function() {
  skip_on_cran()
  tryCatch(
    install_needed("dust2", "compile", FALSE),
    error = function(e) testthat::skip(conditionMessage(e)))
}


logistic_analytic <- function(r, k, times, y0) {
  len <- max(length(r), length(k), length(y0))
  y <- vapply(times, function(t) k / (1 + (k / y0 - 1) * exp(-r * t)),
              numeric(len))
  rbind(y, colSums(y))
}


local_sir_generator <- function() {
  skip_for_compilation()
  gen <- getOption("dust.testing.local_sir_generator", NULL)
  if (!is.null(gen)) {
    return(gen)
  }
  code <- gsub("sir", "mysir", readLines(dust2_file("examples/sir.cpp")))
  filename <- tempfile(fileext = ".cpp")
  writeLines(code, filename)
  gen <- dust_compile(filename, quiet = TRUE, debug = TRUE)
  options(dust.testing.local_sir_generator = gen)
  gen
}


example_bounce <- function(t) {
  skip_if_not_installed("deSolve")
  ball <- function(t, y, parms) {
    dy1 <- y[2]
    dy2 <- -9.8
    list(c(dy1, dy2))
  }
  yini <- c(height = 0, velocity = 10)
  rootfunc <- function(t, y, parms) {
    y[1]
  }
  eventfunc <- function(t, y, parms) {
    y[1] <- 0
    y[2] <- -0.9 * y[2]
    y
  }
  deSolve::ode(times = t, y = yini, func = ball,
               parms = NULL, rootfunc = rootfunc,
               events = list(func = eventfunc, root = TRUE))
}
