test_that("can run system with roots and events", {
  gen <- dust_compile("examples/event.cpp", quiet = TRUE, debug = TRUE)

  sys <- dust_system_create(gen)
  dust_system_set_state_initial(sys)
  t <- seq(0, 6, length.out = 500)
  y <- dust_system_simulate(sys, t)

  ## This is realy not great, but at this point I don't know who is
  ## worse.  Qualitatively we're about right though and I need to go
  ## through and compare with an analytic solution as here we've got
  ## two different approximations and a nonlinear system that is
  ## acumulating error.
  cmp <- example_bounce(t)
  expect_equal(y[1, ], cmp[, 2], tolerance = 1e-2)
})
