test_that("can run simple sir model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust2_cpu_sir_alloc(pars, 0, 1, 10, 0, 42, FALSE)

  expect_type(obj[[1]], "externalptr")
  expect_equal(obj[[2]], 5)

  ptr <- obj[[1]]
  expect_type(dust2_cpu_sir_rng_state(ptr), "raw")
  expect_length(dust2_cpu_sir_rng_state(ptr), 32 * 10)

  expect_equal(dust2_cpu_sir_state(ptr), rep(0, 10 * 5))
  expect_equal(dust2_cpu_sir_time(ptr), 0)

  expect_null(dust2_cpu_sir_set_state_initial(ptr))
  s0 <- dust2_cpu_sir_state(ptr)
  expect_equal(s0, rep(c(990, 10, 0, 0, 0), 10))

  expect_null(dust2_cpu_sir_run_steps(ptr, 30))
  s1 <- matrix(dust2_cpu_sir_state(ptr), 5)
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
    dust2_cpu_sir_compare_data(ptr, d),
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
    dust2_cpu_sir_compare_data(ptr, d),
    rep(0, 10))
  expect_equal(dust2_cpu_sir_rng_state(ptr), r$state())
})
