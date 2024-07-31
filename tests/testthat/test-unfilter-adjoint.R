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


test_that("can't compute adjoint where it was not enabled in the unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))

  ## TODO: these errors end up quite different, but that's
  ## unavoidable?
  obj <- dust_unfilter_create(sir(), time_start, time, data)
  expect_error(dust_unfilter_last_gradient(obj),
               "Gradient is not current")
  ll1 <- dust_unfilter_run(obj, pars)
  expect_error(dust_unfilter_last_gradient(obj),
               "System was not run with 'adjoint = TRUE'")
})


test_that("can compute multiple gradients at once", {
  base <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  n_groups <- 4
  pars <- lapply(seq_len(n_groups),
                 function(i) modifyList(base, list(beta = i * 0.1)))

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) {
    lapply(seq_len(n_groups), function(j) list(incidence = 2 * (i - 1) + j))
  })
  dt <- 1

  obj1 <- dust_unfilter_create(sir(), time_start, time, data, n_groups = 4)
  obj2 <- lapply(1:4, function(i) {
    data_i <- lapply(data, function(x) x[[i]])
    dust_unfilter_create(sir(), time_start, time, data_i)
  })

  ll1 <- dust_unfilter_run(obj1, pars)
  ll2 <- dust_unfilter_run(obj1, pars, adjoint = TRUE)
  expect_equal(ll1, ll2)

  cmp <- vapply(1:4, function(i) {
    dust_unfilter_run(obj2[[i]], pars[[i]], adjoint = TRUE)
    dust_unfilter_last_gradient(obj2[[i]])
  }, numeric(3))

  expect_equal(dust_unfilter_last_gradient(obj1), cmp)
})


test_that("can save history while running unfilter with adjoint", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))

  obj <- dust_unfilter_create(sir(), time_start, time, data)
  ll1 <- dust_unfilter_run(obj, pars, save_history = TRUE)
  h1 <- dust_unfilter_last_history(obj)

  ll2 <- dust_unfilter_run(obj, pars, adjoint = TRUE, save_history = TRUE)
  h2 <- dust_unfilter_last_history(obj)
  expect_identical(ll2, ll1)
  expect_identical(h2, h1)
})


test_that("can run the adjoint model with an odd number of steps", {
  time_start <- 0
  time <- c(3, 6, 9, 12)
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
