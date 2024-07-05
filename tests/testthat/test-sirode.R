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
