test_that("can run an unfilter", {
  base <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  pars1 <- list(beta = 0.1, gamma = 0.2, I0 = 10)
  pars2 <- list(beta = 0.2, gamma = 0.2, I0 = 10)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1

  ## Manually compute likelihood:
  f <- function(pars) {
    base[names(pars)] <- pars
    obj <- dust2_cpu_sir_alloc(base, time_start, dt, 1, 0, NULL, TRUE)
    ptr <- obj[[1]]
    dust2_cpu_sir_set_state_initial(ptr)
    incidence <- numeric(length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust2_cpu_sir_run_steps(ptr, round((time[i] - time0[i]) / dt))
      incidence[i] <- dust2_cpu_sir_state(ptr, FALSE)[5, , drop = TRUE]
    }
    sum(dpois(1:4, incidence + 1e-6, log = TRUE))
  }

  obj <- dust2_cpu_sir_unfilter_alloc(base, time_start, time, dt, data, 0)
  ptr <- obj[[1]]
  expect_equal(dust2_cpu_sir_unfilter_run(ptr, NULL, FALSE), f(pars1))

  expect_equal(dust2_cpu_sir_unfilter_run(ptr, pars1, FALSE), f(pars1))
  expect_equal(dust2_cpu_sir_unfilter_run(ptr, pars2, FALSE), f(pars2))
})
