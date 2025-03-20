test_that("can run system with roots and events", {
  gen <- dust_compile("examples/event.cpp", quiet = TRUE, debug = TRUE)

  control <- dust_ode_control(debug_record_step_times = TRUE, save_history = TRUE)
  sys <- dust_system_create(gen, ode_control = control)
  dust_system_set_state_initial(sys)

  ## Use relatively few points for output here as this exacerbates
  ## problems, even though the solution looks silly.
  t <- seq(0, 6, length.out = 60)
  y <- dust_system_simulate(sys, t)
  cmp <- example_bounce_analytic(t)

  info <- dust_system_internals(sys, include_history = TRUE)

  ## Find all roots:
  r <- info$events[[1]]
  expect_equal(nrow(r), 3)
  expect_equal(r$time, cmp$roots, tolerance = 1e-6)
  expect_equal(r$index, rep(1L, 3))
  expect_equal(r$sign, rep(-1, 3))

  ## Stop at all roots:
  h <- info$history[[1]]
  expect_true(all(r$time %in% h$t0))
  expect_true(all(r$time %in% h$t1))

  ## Overall solution:
  expect_equal(y[1, ], cmp$y, tolerance = 1e-6)
})


test_that("can run events with events in time", {
  gen <- dust_compile("examples/event-time.cpp", quiet = FALSE, debug = TRUE)

  # A mix of times that we will hit exactly and bracket
  pars <- list(r = 0.2, n = 3, t_change = c(2, 5.1234, 7), delta = rnorm(3))
  control <- dust_ode_control(
    debug_record_step_times = TRUE,
    save_history = TRUE
  )
  sys <- dust_system_create(gen, pars, ode_control = control)

  t <- seq(0, 10, length.out = 101)
  y <- dust_system_simulate(sys, t)

  info <- dust_system_internals(sys, include_history = TRUE)
  expect_equal(
    info$events[[1]],
    data_frame(time = pars$t_change, index = 1:3, sign = 1)
  )
  expect_equal(
    drop(y),
    t * 0.2 + c(0, cumsum(pars$delta))[findInterval(t, pars$t_change) + 1])
})


test_that("can cope with coincident events", {
  gen <- dust_compile("examples/event-time.cpp", quiet = FALSE, debug = TRUE)

  pars <- list(r = 0.2, n = 3, t_change = c(2, 2, 3), delta = c(1, 3, 5))
  control <- dust_ode_control(
    debug_record_step_times = TRUE,
    save_history = TRUE
  )
  sys <- dust_system_create(gen, pars, ode_control = control)

  t <- seq(0, 10, length.out = 101)
  y <- dust_system_simulate(sys, t)

  info <- dust_system_internals(sys, include_history = TRUE)
  expect_equal(
    info$events[[1]],
    data_frame(time = pars$t_change, index = 1:3, sign = 1)
  )
  expect_equal(
    drop(y),
    t * 0.2 + c(0, cumsum(pars$delta))[findInterval(t, pars$t_change) + 1])
})
