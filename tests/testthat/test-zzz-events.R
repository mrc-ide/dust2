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
  expect_equal(r$index, rep(0L, 3))
  expect_equal(r$sign, rep(-1, 3))

  ## Stop at all roots:
  h <- info$history[[1]]
  expect_true(all(r$time %in% h$t0))
  expect_true(all(r$time %in% h$t1))

  ## Overall solution:
  expect_equal(y[1, ], cmp$y, tolerance = 1e-6)
})
