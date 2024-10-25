test_that("can run simple logistic system", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  obj <- dust_system_create(logistic(), pars, n_particles = 1,
                            preserve_particle_dimension = TRUE,
                            deterministic = TRUE)
  expect_s3_class(obj, "dust_system")

  expect_equal(dust_system_state(obj), matrix(0, 4, 1))
  expect_equal(dust_system_time(obj), 0)

  dust_system_set_state_initial(obj)
  expect_equal(dust_system_state(obj), matrix(c(1, 1, 1, 3), 4, 1))

  dust_system_run_to_time(obj, 10)
  s <- dust_system_state(obj)
  expect_equal(
    s,
    logistic_analytic(pars$r, pars$K, 10, rep(1, 3)),
    tolerance = 1e-6)

  expect_identical(dim(obj), c(4L, 1L))
  expect_equal(obj$time_control$ode_control, dust_ode_control())
})


test_that("can extract statistics from solver", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  obj <- dust_system_create(logistic(), pars, n_particles = 1,
                            preserve_particle_dimension = TRUE,
                            deterministic = TRUE)

  d <- dust_system_internals(obj)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 1)
  expect_equal(d$particle, 1)
  expect_equal(d$dydt, I(list(c(0, 0, 0))))
  expect_equal(d$step_times, I(list(numeric())))
  expect_equal(d$step_size, 0)
  expect_equal(d$error, 0)
  expect_equal(d$n_steps, 0)
  expect_equal(d$n_steps_accepted, 0)
  expect_equal(d$n_steps_rejected, 0)

  dust_system_set_state_initial(obj)
  dust_system_run_to_time(obj, 10)

  d <- dust_system_internals(obj)
  expect_gt(d$n_steps, 0)
  expect_gt(d$n_steps_accepted, 0)
  ## Step times off by default:
  expect_equal(d$step_times, I(list(numeric())))
})


test_that("can set system state from a vector", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  obj <- dust_system_create(logistic(), pars, n_particles = 1,
                            preserve_particle_dimension = TRUE,
                            deterministic = TRUE)
  expect_s3_class(obj, "dust_system")

  y0 <- matrix(runif(3, max = 20), 3, 1)
  y0 <- rbind(y0, colSums(y0))
  dust_system_set_state(obj, y0)
  expect_equal(dust_system_state(obj), y0)

  dust_system_run_to_time(obj, 10)
  s <- dust_system_state(obj)
  expect_equal(
    s,
    logistic_analytic(pars$r, pars$K, 10, y0[1:3, , drop = FALSE]),
    tolerance = 1e-6)
})


test_that("require an integer value", {
  expect_error(
    dust_system_create(logistic(), list(r = 1, K = 1), n_particles = 1),
    "A value is expected for 'n'")
})


test_that("can set time", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  obj <- dust_system_create(logistic(), pars, n_particles = 10)
  expect_equal(dust_system_time(obj), 0)
  expect_null(dust_system_set_time(obj, 4))
  expect_equal(dust_system_time(obj), 4)
  expect_null(dust_system_set_time(obj, 0))
  expect_equal(dust_system_time(obj), 0)
  expect_null(dust_system_set_time(obj, 0.5))
  expect_equal(dust_system_time(obj), 0.5)
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
  obj <- dust_system_create(logistic(), pars1, n_particles = 1, seed = 42,
                            preserve_particle_dimension = TRUE)

  dust_system_set_state_initial(obj)
  expect_null(dust_system_run_to_time(obj, 5))
  s1 <- dust_system_state(obj)

  expect_null(dust_system_update_pars(obj, pars2))
  expect_null(dust_system_run_to_time(obj, 8))
  s2 <- dust_system_state(obj)

  expect_equal(s1, logistic_analytic(pars1$r, pars1$K, 5, rep(1, 3)),
               tolerance = 1e-6)
  expect_equal(s2, logistic_analytic(pars2$r, pars1$K, 3, s1[1:3, ]),
               tolerance = 1e-4)
})


test_that("can accept a vector of parameters", {
  expect_error(
    dust_system_create(logistic(), list(n = 1, r = 1), n_particles = 1),
    "A value is expected for 'K'")
  expect_error(
    dust_system_create(logistic(), list(n = 1, r = 1, K = c(1, 2)),
                       n_particles = 1),
    "Expected 'K' to have length 1, but it had length 2")
  expect_error(
    dust_system_create(logistic(), list(n = 1, r = 1, K = "a"),
                       n_particles = 1),
    "'K' must be numeric")
})


test_that("can convert integer vectors to numeric", {
  pars1 <- list(n = 3, r = c(0.1, 0.2, 0.3), K = c(100L, 200L, 300L))
  pars2 <- list(n = 3, r = c(0.1, 0.2, 0.3), K = c(100, 200, 300))
  obj1 <- dust_system_create(logistic(), pars1, n_particles = 1,
                             preserve_particle_dimension = TRUE)
  obj2 <- dust_system_create(logistic(), pars2, n_particles = 1,
                             preserve_particle_dimension = TRUE)
  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)
  dust_system_run_to_time(obj1, 10)
  dust_system_run_to_time(obj2, 10)
  expect_identical(dust_system_state(obj1), dust_system_state(obj2))
})


test_that("can reorder particles", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  s <- matrix(runif(30), 3, 10)
  s <- rbind(s, colSums(s))

  obj1 <- dust_system_create(logistic(), pars, n_particles = 10,
                             deterministic = TRUE)
  obj2 <- dust_system_create(logistic(), pars, n_particles = 10,
                             deterministic = TRUE)
  dust_system_set_state(obj1, s)
  dust_system_set_state(obj2, s)

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
  s <- rbind(s, colSums(s))
  dust_system_set_state(obj1, s)
  dust_system_set_state(obj2, s)

  t <- 1:4
  y <- dust_system_simulate(obj1, t)
  expect_equal(dim(y), c(4, 10, 4))

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
    matrix(logistic_analytic(pars[[1]]$r, pars[[1]]$K, 10, rep(1, 3)), 4, 10),
    tolerance = 1e-6)
  expect_equal(
    s[, , 2],
    matrix(logistic_analytic(pars[[2]]$r, pars[[2]]$K, 10, rep(1, 3)), 4, 10),
    tolerance = 1e-6)
})


test_that("can set state into grouped system", {
  pars <- list(
    list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3)),
    list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(200, 3)))

  obj <- dust_system_create(logistic(), pars, n_particles = 10, n_groups = 2)

  s <- array(runif(3 * 10 * 2), c(3, 10, 2))
  dust_system_set_state(obj, s, index_state = 1:3)

  expect_equal(dust_system_state(obj, index_state = 1:3), s)
  expect_equal(dust_system_state(obj, index_state = 4),
               array(apply(s, 2:3, sum), c(1, 10, 2)))

  s <- array(runif(4 * 2), c(4, 1, 2))
  dust_system_set_state(obj, s)
  expect_equal(dust_system_state(obj), s[, rep(1, 10), ])

  s <- array(runif(4 * 10), c(4, 10, 1))
  dust_system_set_state(obj, s)
  expect_equal(dust_system_state(obj), s[, , c(1, 1)])

  s <- array(runif(4), c(4, 1, 1))
  dust_system_set_state(obj, s)
  expect_equal(dust_system_state(obj), s[, rep(1, 10), c(1, 1)])
})


test_that("reject dt as an argument when creating logistic model", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  expect_error(
    dust_system_create(logistic(), pars, dt = 1, n_particles = 1),
    "Can't use 'dt' with continuous-time systems")
})


test_that("can set ode control", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  ctl1 <- dust_ode_control(atol = 1e-8, rtol = 1e-8)
  ctl2 <- dust_ode_control(atol = 1e-3, rtol = 1e-3)
  obj1 <- dust_system_create(logistic(), pars, n_particles = 1,
                             preserve_particle_dimension = TRUE,
                             ode_control = ctl1)
  obj2 <- dust_system_create(logistic(), pars, n_particles = 1,
                             preserve_particle_dimension = TRUE,
                             ode_control = ctl2)
  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)
  dust_system_run_to_time(obj1, 10)
  dust_system_run_to_time(obj2, 10)
  s1 <- dust_system_state(obj1)
  s2 <- dust_system_state(obj2)
  expect_false(identical(s1, s2))

  cmp <- logistic_analytic(pars$r, pars$K, 10, rep(1, 3))
  expect_equal(s1, cmp, tolerance = 1e-7)
  expect_equal(s2, cmp, tolerance = 1e-3)
  expect_true(all(abs(s1 - cmp) < abs(s2 - cmp)))

  expect_equal(obj1$time_control$ode_control, ctl1)
  expect_equal(obj2$time_control$ode_control, ctl2)
})


test_that("can error if too many steps taken", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  ctl <- dust_ode_control(max_steps = 2)
  obj <- dust_system_create(logistic(), pars, n_particles = 1,
                            preserve_particle_dimension = TRUE,
                            ode_control = ctl)
  dust_system_set_state_initial(obj)
  expect_error(dust_system_run_to_time(obj, 10),
               "too many steps")
})


test_that("can error if steps are too small", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  ctl <- dust_ode_control(step_size_min = 10)
  obj <- dust_system_create(logistic(), pars, n_particles = 1,
                            preserve_particle_dimension = TRUE,
                            ode_control = ctl)
  dust_system_set_state_initial(obj)
  expect_error(dust_system_run_to_time(obj, 10),
               "step too small")
})


test_that("can error if initial step size calculation fails", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  obj <- dust_system_create(logistic(), pars, n_particles = 1,
                            preserve_particle_dimension = TRUE,
                            deterministic = TRUE)
  y0 <- matrix(NA_real_, 4, 1)
  expect_error(dust_system_set_state(obj, y0),
               "Initial step size was not finite")
})


test_that("can save step times for debugging", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  ctl <- dust_ode_control(debug_record_step_times = TRUE)
  obj <- dust_system_create(logistic(), pars, n_particles = 1,
                            preserve_particle_dimension = TRUE,
                            deterministic = TRUE, ode_control = ctl)

  dust_system_set_state_initial(obj)
  dust_system_run_to_time(obj, 10)

  d <- dust_system_internals(obj)
  expect_type(d$step_times, "list")
  expect_length(d$step_times[[1]], d$n_steps + 1)
  expect_equal(d$step_times[[1]][[1]], 0)
  expect_equal(d$step_times[[1]][[d$n_steps + 1]], 10)
  expect_false(any(diff(d$step_times[[1]]) <= 0))
})




test_that("can save history", {
  pars <- list(n = 3, r = c(0.1, 0.2, 0.3), K = rep(100, 3))
  ctl <- dust_ode_control(save_history = TRUE, debug_record_step_times = TRUE)
  sys <- dust_system_create(logistic(), pars, n_particles = 1,
                            preserve_particle_dimension = TRUE,
                            deterministic = TRUE, ode_control = ctl)
  dust_system_set_state_initial(sys)
  dust_system_run_to_time(sys, 10)
  d <- dust_system_internals(sys,
                             include_coefficients = TRUE,
                             include_history = TRUE)
  expect_true("history" %in% names(d))
  history <- d$history[[1]]
  n_steps <- d$n_steps_accepted
  expect_length(history$time, n_steps)
  expect_equal(history$time, d$step_times[[1]][-(n_steps + 1)])
  expect_equal(history$size, diff(d$step_times[[1]]))
  expect_equal(history$coefficients[, , n_steps], d$coefficients[[1]])
})
