test_that("can run sirode model", {
  pars <- list(beta = 0.2, gamma = 0.1, N = 1000, I0 = 10, exp_noise = 1e6)
  obj <- dust_system_create(sirode(), pars, n_particles = 1,
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
                               ode_control = ctl, deterministic = TRUE)
    dust_system_set_state_initial(obj1)
    dust_system_run_to_time(obj1, t)
    list(times = dust_system_internals(obj1)$step_times[[1]],
         state = dust_system_state(obj1))
  }

  res <- lapply(1:10, f)

  obj2 <- dust_system_create(sirode(), pars, n_particles = 1,
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
