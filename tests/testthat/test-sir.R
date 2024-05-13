test_that("can run simple sir model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust2_cpu_sir_alloc(pars, 0, 1, 10, 0, 42, FALSE)

  expect_type(obj[[1]], "externalptr")
  expect_equal(obj[[2]], 5)

  ptr <- obj[[1]]
  expect_type(dust2_cpu_sir_rng_state(ptr), "raw")
  expect_length(dust2_cpu_sir_rng_state(ptr), 32 * 10)

  expect_equal(dust2_cpu_sir_state(ptr, FALSE), matrix(0, 5, 10))
  expect_equal(dust2_cpu_sir_time(ptr), 0)

  expect_null(dust2_cpu_sir_set_state_initial(ptr))
  s0 <- dust2_cpu_sir_state(ptr, FALSE)
  expect_equal(s0, matrix(c(990, 10, 0, 0, 0), 5, 10))

  expect_null(dust2_cpu_sir_run_steps(ptr, 30))
  s1 <- dust2_cpu_sir_state(ptr, FALSE)
  expect_true(all(s1[1, ] < 990))
  expect_true(all(s1[3, ] > 0))
  expect_true(all(s1[4, ] > 0))
})


test_that("can compare to data", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 0.5)
  obj <- dust2_cpu_sir_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]

  s <- rbind(0, 0, 0, 0, rpois(10, 30))
  dust2_cpu_sir_set_state(ptr, s)
  d <- list(incidence = 30)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  eps <- drop(r$exponential(1, 0.5))

  expect_equal(
    dust2_cpu_sir_compare_data(ptr, d, FALSE),
    dpois(30, s[5, ] + eps, log = TRUE))
})


test_that("can compare to data when missing", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 0.5)
  obj <- dust2_cpu_sir_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]

  s <- rbind(0, 0, 0, 0, rpois(10, 30))
  dust2_cpu_sir_set_state(ptr, s)
  d <- list(incidence = NA_real_)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  expect_equal(
    dust2_cpu_sir_compare_data(ptr, d, FALSE),
    rep(0, 10))
  expect_equal(dust2_cpu_sir_rng_state(ptr), r$state())
})


test_that("can compare against multple parameter groups at once", {
  pars <- lapply(1:4, function(i) {
    list(beta = 0.1 * i, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 10^i)
  })
  obj <- dust2_cpu_sir_alloc(pars, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]

  s <- dust2_cpu_sir_state(ptr, TRUE)
  s[5, , ] <- rpois(10, 30)
  dust2_cpu_sir_set_state(ptr, s)

  d <- lapply(1:4, function(i) list(incidence = 30 + i))
  res <- dust2_cpu_sir_compare_data(ptr, d, TRUE)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10 * 4)
  rate <- rep(10^(1:4), each = 10)
  eps <- matrix(r$exponential(1, 1) /  rate, 10, 4)
  expect_equal(
    res,
    matrix(dpois(rep(31:34, each = 10), s[5, , ] + eps, log = TRUE), 10, 4))

  expect_error(
    dust2_cpu_sir_compare_data(ptr, d, FALSE),
    "Can't compare with grouped = FALSE with more than one group")
})


test_that("validate data size on compare", {
  pars <- lapply(1:4, function(i) {
    list(beta = 0.1 * i, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 10^i)
  })
  obj <- dust2_cpu_sir_alloc(pars, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]
  expect_error(
    dust2_cpu_sir_compare_data(ptr, vector("list", 3), TRUE),
    "'data' must have length 4")
})


test_that("can update parameters", {
  base <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  update <- list(beta = 0.3, gamma = 0.2)
  combined <- modifyList(base, update)

  obj1 <- dust2_cpu_sir_alloc(base, 0, 1, 10, 0, 42, FALSE)
  ptr1 <- obj1[[1]]
  expect_null(dust2_cpu_sir_update_pars(ptr1, update, FALSE))
  expect_null(dust2_cpu_sir_run_steps(ptr1, 10))

  obj2 <- dust2_cpu_sir_alloc(combined, 0, 1, 10, 0, 42, FALSE)
  ptr2 <- obj2[[1]]
  expect_null(dust2_cpu_sir_run_steps(ptr2, 10))

  expect_equal(
    dust2_cpu_sir_state(ptr2, FALSE),
    dust2_cpu_sir_state(ptr1, FALSE))
})
