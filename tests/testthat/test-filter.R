test_that("can create an unfilter object", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  obj <- dust2_cpu_sir_unfilter_alloc(pars, time_start, time, dt, data, 0)
  ptr <- obj[[1]]
  ll <- dust2_cpu_sir_unfilter_run(ptr, NULL)
})
