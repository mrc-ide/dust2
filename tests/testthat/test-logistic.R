test_that("can run simple logistic system", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  obj <- dust_system_create(logistic(), pars, n_particles = 1,
                            deterministic = TRUE)
  expect_s3_class(obj, "dust_system")

  expect_equal(dust_system_state(obj), matrix(0, 3, 1))
  expect_equal(dust_system_time(obj), 0)

  dust_system_set_state_initial(obj)
  expect_equal(dust_system_state(obj), matrix(1, 3, 1))

  dust_system_run_to_time(obj, 10)
  s <- dust_system_state(obj)
  expect_equal(
    s,
    logistic_analytic(pars$r, pars$K, 10, rep(1, 3)),
    tolerance = 1e-6)

  expect_identical(dim(obj), c(3L, 1L))
})


test_that("can set system state from a vector", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  obj <- dust_system_create(logistic(), pars, n_particles = 1,
                            deterministic = TRUE)
  expect_s3_class(obj, "dust_system")

  y0 <- matrix(runif(3, max = 20), 3, 1)
  dust_system_set_state(obj, y0)
  expect_equal(dust_system_state(obj), y0)

  dust_system_run_to_time(obj, 10)
  s <- dust_system_state(obj)
  expect_equal(
    s,
    logistic_analytic(pars$r, pars$K, 10, y0),
    tolerance = 1e-6)
})


test_that("can set time", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  obj <- dust_system_create(logistic(), pars, n_particles = 10)
  expect_equal(dust_system_time(obj), 0)
  expect_null(dust_system_set_time(obj, 4))
  expect_equal(dust_system_time(obj), 4)
  expect_null(dust_system_set_time(obj, 0))
  expect_equal(dust_system_time(obj), 0)
  expect_error(dust_system_set_time(obj, 0.5),
               "Expected 'time' to be integer-like")
})


test_that("can set rng state", {
  pars <- list(n = 1, r = 1, K = 100)

  obj1 <- dust_system_create(logistic(), pars, n_particles = 10, seed = 42)
  obj2 <- dust_system_create(logistic(), pars, n_particles = 10, seed = 43)

  expect_false(identical(dust_system_rng_state(obj1),
                         dust_system_rng_state(obj2)))

  expect_null(dust_system_set_rng_state(obj2, dust_system_rng_state(obj1)))
  expect_identical(dust_system_rng_state(obj1), dust_system_rng_state(obj2))
  expect_identical(dust_system_state(obj1), dust_system_state(obj2))
})


test_that("can update parameters", {
  pars1 <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  pars2 <- list(r = pars1$r * 5)
  obj <- dust_system_create(logistic(), pars1, n_particles = 1, seed = 42)

  dust_system_set_state_initial(obj)
  expect_null(dust_system_run_to_time(obj, 5))
  s1 <- dust_system_state(obj)

  expect_null(dust_system_update_pars(obj, pars2))
  expect_null(dust_system_run_to_time(obj, 8))
  s2 <- dust_system_state(obj)

  expect_equal(s1, logistic_analytic(pars1$r, pars1$K, 5, rep(1, 3)),
               tolerance = 1e-6)
  expect_equal(s2, logistic_analytic(pars2$r, pars1$K, 3, s1),
               tolerance = 1e-4)
})
