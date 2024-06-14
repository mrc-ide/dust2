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
  expect_equal(obj$ode_control, dust_ode_control())
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


test_that("can accept a vector of parameters", {
  expect_error(
    dust_system_create(logistic(), list(n = 1, r = 1), n_particles = 1),
    "A value is expected for 'K'")
  expect_error(
    dust_system_create(logistic(), list(n = 1, r = 1, K = c(1, 2)),
                       n_particles = 1),
    "'K' must have length 1")
  expect_error(
    dust_system_create(logistic(), list(n = 1, r = 1, K = "a"),
                       n_particles = 1),
    "'K' must be numeric")
})


test_that("can convert integer vectors to numeric", {
  pars1 <- list(n = 3, r = c(0.1, 0.2, 0.3), K = c(100L, 200L, 300L))
  pars2 <- list(n = 3, r = c(0.1, 0.2, 0.3), K = c(100, 200, 300))
  obj1 <- dust_system_create(logistic(), pars1, n_particles = 1)
  obj2 <- dust_system_create(logistic(), pars2, n_particles = 1)
  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)
  dust_system_run_to_time(obj1, 10)
  dust_system_run_to_time(obj2, 10)
  expect_identical(dust_system_state(obj1), dust_system_state(obj2))
})


test_that("can reorder particles", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  s <- matrix(runif(30), 3, 10)

  obj1 <- dust_system_create(logistic(), pars, n_particles = 10,
                             deterministic = TRUE)
  obj2 <- dust_system_create(logistic(), pars, n_particles = 10,
                             deterministic = TRUE)
  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)

  dust_system_run_to_time(obj1, 3)
  dust_system_run_to_time(obj1, 10)

  dust_system_run_to_time(obj2, 3)
  dust_system_reorder(obj2, 10:1)
  dust_system_run_to_time(obj2, 10)

  s1 <- dust_system_state(obj1)
  s2 <- dust_system_state(obj2)
  expect_identical(s2, s1[, 10:1])
})


## Later, we might change simulate to take advantage of dense output,
## but for now we stop at every time exactly.
test_that("can simulate", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  obj1 <- dust_system_create(logistic(), pars, n_particles = 10,
                             deterministic = TRUE)
  obj2 <- dust_system_create(logistic(), pars, n_particles = 10,
                             deterministic = TRUE)
  s <- matrix(runif(30), 3, 10)
  dust_system_set_state(obj1, s)
  dust_system_set_state(obj2, s)

  t <- 1:4
  y <- dust_system_simulate(obj1, t)
  expect_equal(dim(y), c(3, 10, 4))

  dust_system_run_to_time(obj2, 1)
  expect_identical(dust_system_state(obj2), y[, , 1])
  dust_system_run_to_time(obj2, 2)
  expect_identical(dust_system_state(obj2), y[, , 2])
})


test_that("can create grouped system", {
  pars <- list(
    list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3)),
    list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(200, 3)))

  obj <- dust_system_create(logistic(), pars, n_particles = 10, n_groups = 2)
  dust_system_set_state_initial(obj)
  dust_system_run_to_time(obj, 10)
  s <- dust_system_state(obj)

  expect_equal(
    s[, , 1],
    matrix(logistic_analytic(pars[[1]]$r, pars[[1]]$K, 10, rep(1, 3)), 3, 10),
    tolerance = 1e-6)
  expect_equal(
    s[, , 2],
    matrix(logistic_analytic(pars[[2]]$r, pars[[2]]$K, 10, rep(1, 3)), 3, 10),
    tolerance = 1e-6)
})


test_that("can set state into grouped system", {
  pars <- list(
    list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3)),
    list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(200, 3)))

  obj <- dust_system_create(logistic(), pars, n_particles = 10, n_groups = 2)

  s <- array(runif(3 * 10 * 2), c(3, 10, 2))
  dust_system_set_state(obj, s)
  expect_equal(dust_system_state(obj), s)

  s <- array(runif(3 * 2), c(3, 1, 2))
  dust_system_set_state(obj, s)
  expect_equal(dust_system_state(obj), s[, rep(1, 10), ])

  s <- array(runif(3 * 10), c(3, 10, 1))
  dust_system_set_state(obj, s)
  expect_equal(dust_system_state(obj), s[, , c(1, 1)])

  s <- array(runif(3), c(3, 1, 1))
  dust_system_set_state(obj, s)
  expect_equal(dust_system_state(obj), s[, rep(1, 10), c(1, 1)])
})


test_that("reject dt as an argument when creating logistic model", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  expect_error(
    dust_system_create(logistic(), pars, dt = 1, n_particles = 1),
    "Can't use 'dt' with continuous-time systems")
})
