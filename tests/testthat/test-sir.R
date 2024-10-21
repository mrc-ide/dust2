test_that("can run simple sir system", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)

  expect_type(obj$ptr, "externalptr")
  expect_equal(obj$n_state, 5)

  expect_equal(obj$packing_state,
               list(S = integer(), I = integer(), R = integer(),
                    cases_cumul = integer(), cases_inc = integer()))
  expect_s3_class(obj$packer_state, "monty_packer")

  expect_type(dust_system_rng_state(obj), "raw")
  expect_length(dust_system_rng_state(obj), 32 * 10)

  expect_equal(dust_system_state(obj), matrix(0, 5, 10))
  expect_equal(dust_system_time(obj), 0)

  expect_null(dust_system_set_state_initial(obj))
  s0 <- dust_system_state(obj)
  expect_equal(s0, matrix(c(990, 10, 0, 0, 0), 5, 10))

  expect_null(dust_system_run_to_time(obj, 30))
  s1 <- dust_system_state(obj)
  expect_true(all(s1[1, ] < 990))
  expect_true(all(s1[3, ] > 0))
  expect_true(all(s1[4, ] > 0))
})


test_that("can compare to data", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 0.5)
  obj <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)

  s <- rbind(0, 0, 0, 0, rpois(10, 30))
  dust_system_set_state(obj, s)
  d <- list(incidence = 30)

  r <- monty::monty_rng$new(seed = 42, n_streams = 10)
  eps <- drop(r$exponential_rate(1, 0.5))

  expect_equal(
    dust_system_compare_data(obj, d),
    dpois(30, s[5, ] + eps, log = TRUE))
})


test_that("can compare to data when missing", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 0.5)
  obj <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)

  s <- rbind(0, 0, 0, 0, rpois(10, 30))
  dust_system_set_state(obj, s)
  d <- list(incidence = NA_real_)

  r <- monty::monty_rng$new(seed = 42, n_streams = 10)
  expect_equal(
    dust_system_compare_data(obj, d),
    rep(0, 10))
  expect_equal(dust_system_rng_state(obj), r$state())
})


test_that("can compare against multple parameter groups at once", {
  pars <- lapply(1:4, function(i) {
    list(beta = 0.1 * i, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 10^i)
  })
  obj <- dust_system_create(sir(), pars, n_particles = 10, n_groups = 4,
                           seed = 42)

  s <- dust_system_state(obj)
  s[5, , ] <- rpois(10, 30)
  dust_system_set_state(obj, s)

  d <- lapply(1:4, function(i) list(incidence = 30 + i))
  res <- dust_system_compare_data(obj, d)

  r <- monty::monty_rng$new(seed = 42, n_streams = 10 * 4)
  rate <- rep(10^(1:4), each = 10)
  eps <- matrix(r$exponential_rate(1, 1) /  rate, 10, 4)
  expect_equal(
    res,
    matrix(dpois(rep(31:34, each = 10), s[5, , ] + eps, log = TRUE), 10, 4))
})


test_that("validate data size on compare", {
  pars <- lapply(1:4, function(i) {
    list(beta = 0.1 * i, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 10^i)
  })
  obj <- dust_system_create(sir(), pars, n_particles = 10, n_groups = 4)
  expect_error(
    dust_system_compare_data(obj, vector("list", 3)),
    "Expected 'data' to have length 4, but it had length 3")
})


test_that("periodic calculation adds up to cumulative cases", {
  n_time <- 20
  n_particles <- 100

  pars <- list(beta = 1.0, gamma = 0.5, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust_system_create(sir(), pars, n_particles = n_particles,
                            seed = 42)
  dust_system_set_state_initial(obj)
  res <- dust_system_simulate(obj, 0:5)

  ## Cumulative cases never decrease:
  expect_true(all(diff(t(res[4, , ])) >= 0))

  ## Incidence resets somtimes:
  expect_true(any(diff(t(res[5, , ])) < 0))

  expect_equal(apply(res[5, , ], 1, cumsum), t(res[4, , ]))
})


test_that("can reset cases daily", {
  n_time <- 20
  n_particles <- 100

  pars <- list(beta = 1.0, gamma = 0.5, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust_system_create(sir(), pars, n_particles = n_particles,
                           dt = 0.25, seed = 42)
  dust_system_set_state_initial(obj)
  res <- dust_system_simulate(obj, 0:5)

  ## Cumulative cases never decrease:
  expect_true(all(diff(t(res[4, , ])) >= 0))

  ## Incidence resets somtimes:
  expect_true(any(diff(t(res[5, , ])) < 0))

  expect_equal(apply(res[5, , ], 1, cumsum), t(res[4, , ]))
})


test_that("can run sir system to time", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  obj1 <- dust_system_create(sir(), pars, dt = 0.25, n_particles = 10,
                             seed = 42)
  obj2 <- dust_system_create(sir(), pars, dt = 0.25, n_particles = 10,
                             seed = 42)

  dust_system_set_state_initial(obj1)
  expect_null(dust_system_run_to_time(obj1, 10))
  expect_equal(dust_system_time(obj1), 10)
  s1 <- dust_system_state(obj1)

  dust_system_set_state_initial(obj2)
  expect_null(dust_system_run_to_time(obj2, 10))
  expect_equal(dust_system_time(obj2), 10)
  expect_equal(dust_system_state(obj2), s1)
})


test_that("can update parameters", {
  base <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  update <- list(beta = 0.3, gamma = 0.2)
  combined <- modifyList(base, update)

  obj1 <- dust_system_create(sir(), base, n_particles = 10, seed = 42)
  expect_null(dust_system_update_pars(obj1, update))
  expect_null(dust_system_run_to_time(obj1, 10))

  obj2 <- dust_system_create(sir(), combined, n_particles = 10, seed = 42)
  expect_null(dust_system_run_to_time(obj2, 10))

  expect_equal(
    dust_system_state(obj2),
    dust_system_state(obj1))
})


test_that("can run simulation", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  obj1 <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)
  obj2 <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)
  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)

  res <- dust_system_simulate(obj1, 0:20)
  expect_equal(dim(res), c(5, 10, 21))

  expect_equal(res[, , 21], dust_system_state(obj1))
  expect_equal(dust_system_time(obj1), 20)

  expect_equal(res[, , 1], dust_system_state(obj2))
  dust_system_run_to_time(obj2, 10)
  expect_equal(res[, , 11], dust_system_state(obj2))
})


test_that("can run simulation with index", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  obj1 <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)
  obj2 <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)
  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)

  index <- c(2L, 4L)
  res1 <- dust_system_simulate(obj1, 0:20, index)
  res2 <- dust_system_simulate(obj2, 0:20)

  expect_equal(res1, res2[index, , ])
})


test_that("copy names with index", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  obj1 <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)
  obj2 <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)
  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)

  index <- c(I = 2L, cases = 4L)
  res1 <- dust_system_simulate(obj1, 0:20, index)
  res2 <- dust_system_simulate(obj2, 0:20, unname(index))
  expect_equal(dimnames(res1), list(c("I", "cases"), NULL, NULL))

  expect_equal(unname(res1), res2)
})


test_that("can reorder state", {
  obj <- dust_system_create(sir(), list(), n_particles = 10, seed = 42)
  expect_null(dust_system_set_state_initial(obj))
  dust_system_run_to_time(obj, 100)
  s1 <- dust_system_state(obj)
  i <- sample(10, replace = TRUE)
  expect_null(dust_system_reorder(obj, i))
  s2 <- dust_system_state(obj)
  expect_equal(s2, s1[, i, drop = FALSE])
})


test_that("can set rng state", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  obj1 <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)
  obj2 <- dust_system_create(sir(), pars, n_particles = 10, seed = 43)

  expect_false(identical(dust_system_rng_state(obj1),
                         dust_system_rng_state(obj2)))

  expect_null(dust_system_set_rng_state(obj2, dust_system_rng_state(obj1)))
  expect_identical(dust_system_rng_state(obj1), dust_system_rng_state(obj2))
})


test_that("throw nice error on run failure, then prevent access", {
  pars <- list(beta = 0.1, gamma = -0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust_system_create(sir(), pars, n_particles = 10, seed = 42)
  expect_null(dust_system_set_state_initial(obj))

  err <- expect_error(
    expect_null(dust_system_run_to_time(obj, 30)),
    "10 particles reported errors")
  expect_match(conditionMessage(err), "1: Invalid call to binomial with")
  expect_match(conditionMessage(err), "- (and 6 more)", fixed = TRUE)

  expect_error(
    dust_system_run_to_time(obj, 30),
    "Can't currently run this system: errors are pending")
  expect_error(
    dust_system_time(obj),
    "Can't currently get time from this system: errors are pending")
  expect_error(
    dust_system_state(obj),
    "Can't currently get state from this system: errors are pending")
  expect_error(
    dust_system_reorder(obj, 1:10),
    "Can't currently reorder this system: errors are pending")
  expect_error(
    dust_system_set_time(obj, 0),
    "Can't currently set time for this system: errors are pending")
  expect_error(
    dust_system_compare_data(obj, list()),
    "Can't currently compare data for this system: errors are pending")
  expect_error(
    dust_system_simulate(obj, 0:10),
    "Can't currently simulate this system: errors are pending")

  dust_system_update_pars(obj, list(gamma = 0.2))
  dust_system_set_state_initial(obj)
  expect_equal(dust_system_time(obj), 0)
})


test_that("can unpack state", {
  set.seed(1)
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust_system_create(sir(), pars, n_particles = 1)
  dust_system_set_state_initial(obj)
  m <- dust_system_state(obj)
  expect_equal(obj$packer_state$unpack(m),
               list(S = 990, I = 10, R = 0, cases_cumul = 0, cases_inc = 0))
})


test_that("can set state from a vector", {
  set.seed(1)
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  sys <- dust_system_create(sir(), pars, n_particles = 10)
  s <- c(1000, 10, 0, 0, 0)
  dust_system_set_state(sys, s)
  expect_equal(dust_system_state(sys), matrix(s, 5, 10))
})
