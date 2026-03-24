test_that("continuous adjoint gives same likelihood as forward", {
  pars <- list(beta = 0.2, gamma = 0.1, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sirode_adjoint(), time_start, data)
  ll1 <- dust_likelihood_run(obj, pars, adjoint = FALSE)
  ll2 <- dust_likelihood_run(obj, pars, adjoint = TRUE)
  expect_equal(ll2, ll1, tolerance = 1e-6)
})


test_that("continuous adjoint gradient matches finite differences", {
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sirode_adjoint(), time_start, data)
  x <- c(beta = 0.2, gamma = 0.1, I0 = 10)
  ll <- dust_likelihood_run(obj, as.list(x), adjoint = TRUE)
  gr <- dust_likelihood_last_gradient(obj)

  gr_num <- numDeriv::grad(
    function(x) dust_likelihood_run(obj, as.list(x)),
    x, method = "Richardson", method.args = list(r = 6))
  expect_equal(gr_num, gr, tolerance = 1e-4)
})


test_that("continuous adjoint gradient correct at different parameters", {
  time_start <- 0
  data <- data.frame(time = seq(5, 50, by = 5),
                     incidence = c(3, 8, 15, 25, 30, 28, 20, 12, 7, 4))

  obj <- dust_unfilter_create(sirode_adjoint(), time_start, data)

  for (beta in c(0.15, 0.3, 0.5)) {
    x <- c(beta = beta, gamma = 0.1, I0 = 10)
    dust_likelihood_run(obj, as.list(x), adjoint = TRUE)
    gr <- dust_likelihood_last_gradient(obj)

    gr_num <- numDeriv::grad(
      function(x) dust_likelihood_run(obj, as.list(x)),
      x, method = "Richardson", method.args = list(r = 6))
    expect_equal(gr_num, gr, tolerance = 1e-3,
                 label = sprintf("gradient at beta=%g", beta))
  }
})


test_that("continuous adjoint handles NA data correctly", {
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16),
                     incidence = c(1, NA, 3, NA))

  obj <- dust_unfilter_create(sirode_adjoint(), time_start, data)
  x <- c(beta = 0.2, gamma = 0.1, I0 = 10)
  ll <- dust_likelihood_run(obj, as.list(x), adjoint = TRUE)
  gr <- dust_likelihood_last_gradient(obj)

  expect_true(all(is.finite(gr)))

  gr_num <- numDeriv::grad(
    function(x) dust_likelihood_run(obj, as.list(x)),
    x, method = "Richardson", method.args = list(r = 6))
  expect_equal(gr_num, gr, tolerance = 1e-4)
})


test_that("continuous adjoint works via monty interface", {
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sirode_adjoint(), time_start, data)
  packer <- monty::monty_packer(c("beta", "gamma", "I0"),
                                fixed = list(N = 1000, exp_noise = 1e6))
  ll <- dust_likelihood_monty(obj, packer)

  expect_true(ll$properties$has_gradient)

  theta <- c(0.2, 0.1, 10)
  d <- ll$density(theta)
  g <- ll$gradient(theta)

  expect_true(is.finite(d))
  expect_length(g, 3)
  expect_true(all(is.finite(g)))
})


test_that("error if gradient requested without adjoint run", {
  pars <- list(beta = 0.2, gamma = 0.1, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sirode_adjoint(), time_start, data)
  expect_error(dust_likelihood_last_gradient(obj),
               "Gradient is not current")
  dust_likelihood_run(obj, pars, adjoint = FALSE)
  expect_error(dust_likelihood_last_gradient(obj),
               "System was not run with 'adjoint = TRUE'")
})
