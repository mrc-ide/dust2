test_that("can run sirode model", {
  pars <- list(beta = 0.2, gamma = 0.1, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust_system_create(sirode(), pars, n_particles = 1,
                            preserve_particle_dimension = TRUE,
                            ode_control = dust_ode_control(step_size_max = 0.8),
                            deterministic = TRUE)
  dust_system_set_state_initial(obj)
  res <- dust_system_simulate(obj, 0:10)

  ## Cumulative cases never decrease:
  expect_true(all(diff(res[4, , ]) >= 0))

  expect_equal(cumsum(res[5, , ]), res[4, , ])
})


## Rather than the existing loop and double checking
test_that("can cope with multiple resets within a step", {
  ctl <- dust_ode_control(debug_record_step_times = TRUE)
  pars <- list(beta = 0.2, gamma = 0.1, N = 1000, I0 = 10, exp_noise = 1e6)

  f <- function(t) {
    obj1 <- dust_system_create(sirode(), pars, n_particles = 1,
                               preserve_particle_dimension = TRUE,
                               ode_control = ctl, deterministic = TRUE)
    dust_system_set_state_initial(obj1)
    dust_system_run_to_time(obj1, t)
    list(times = dust_system_internals(obj1)$step_times[[1]],
         state = dust_system_state(obj1))
  }

  res <- lapply(1:10, f)

  obj2 <- dust_system_create(sirode(), pars, n_particles = 1,
                             preserve_particle_dimension = TRUE,
                             ode_control = ctl, deterministic = TRUE)
  dust_system_set_state_initial(obj2)
  cmp <- dust_system_simulate(obj2, 1:10)

  ans <- vnapply(res, function(x) x$state[[5]])
  expect_equal(ans, cmp[5, , ], tolerance = 1e-5)

  expect_true(any(vnapply(res, function(x) diff(tail(x$times, 2))) > 1))
})


test_that("can extract interpolation coefficients from internals", {
  pars <- list(beta = 0.2, gamma = 0.1, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust_system_create(sirode(), pars, n_particles = 1,
                            ode_control = dust_ode_control(step_size_max = 0.8),
                            deterministic = TRUE)
  dust_system_set_state_initial(obj)
  res <- dust_system_run_to_time(obj, 10)

  d1 <- dust_system_internals(obj)
  expect_s3_class(d1, "data.frame")
  expect_false("coefficients" %in% names(d1))

  d2 <- dust_system_internals(obj, include_coefficients = TRUE)
  expect_s3_class(d2, "data.frame")
  expect_identical(d2[names(d2) != "coefficients"], d1)
  expect_true("coefficients" %in% names(d2))
  expect_equal(dim(d2$coefficients[[1]]), c(5, 5))
})


test_that("throw nice error on run failure, then prevent access", {
  pars <- list(beta = 0.2, gamma = 0.1, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust_system_create(sirode(), pars, n_particles = 10,
                            ode_control = dust_ode_control(max_steps = 5),
                            deterministic = TRUE)
  dust_system_set_state_initial(obj)

  err <- expect_error(
    dust_system_run_to_time(obj, 30),
    "10 particles reported errors")
  expect_match(conditionMessage(err), "1: too many steps")
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
    dust_system_simulate(obj, 0:10),
    "Can't currently simulate this system: errors are pending")

  dust_system_update_pars(obj, list(gamma = 0.2))
  dust_system_set_state_initial(obj)
  expect_equal(dust_system_time(obj), 0)
})


test_that("can compare to data", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 0.5)
  obj <- dust_system_create(sirode(), pars, 10, deterministic = TRUE)

  s <- rbind(0, 0, 0, 0, rpois(10, 30) + rnorm(10))
  dust_system_set_state(obj, s)
  d <- list(incidence = 30)

  expect_equal(
    dust_system_compare_data(obj, d),
    dpois(30, s[5, ], log = TRUE))
})


test_that("can force solver to stop at times", {
  pars <- list(beta = 0.2, gamma = 0.1, N = 1000, I0 = 10, exp_noise = 1e6)
  ctl1 <- dust_ode_control(critical_times = c(10, 20, 30),
                           debug_record_step_times = TRUE)
  ctl2 <- dust_ode_control(debug_record_step_times = TRUE)

  obj1 <- dust_system_create(sirode, pars, n_particles = 1,
                             ode_control = ctl1, deterministic = TRUE)
  obj2 <- dust_system_create(sirode, pars, n_particles = 1,
                             ode_control = ctl2, deterministic = TRUE)

  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)
  dust_system_run_to_time(obj1, 50)
  dust_system_run_to_time(obj2, 50)

  t1 <- dust_system_internals(obj1)$step_times[[1]]
  t2 <- dust_system_internals(obj2)$step_times[[1]]
  expect_true(all(c(10, 20, 30) %in% t1))
  expect_false(all(c(10, 20, 30) %in% t2))

  expect_equal(dust_system_state(obj1),
               dust_system_state(obj2), tolerance = 1e-6)
})
