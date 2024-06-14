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
