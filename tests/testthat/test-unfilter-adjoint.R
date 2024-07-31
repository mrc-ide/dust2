test_that("can run an unfilter via the adjoint method", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))

  obj <- dust_unfilter_create(sir(), time_start, time, data)
  ll1 <- dust_unfilter_run(obj, pars)
  ll2 <- dust_unfilter_run(obj, pars, adjoint = TRUE)
  expect_identical(ll2, ll1)
})


test_that("can run the adjoint model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))

  obj <- dust_unfilter_create(sir(), time_start, time, data)
  x <- c(beta = 0.1, gamma = 0.2, I0 = 10)
  ll <- dust_unfilter_run(obj, as.list(x), adjoint = TRUE)
  gr <- dust_unfilter_last_gradient(obj)

  ll <- dust_unfilter_run(obj, as.list(x))
  gr_num <- numDeriv::grad(function(x) dust_unfilter_run(obj, as.list(x)), x,
                           method = "Richardson", method.args = list(r = 6))
  expect_equal(gr_num, gr, tolerance = 1e-4)
})
