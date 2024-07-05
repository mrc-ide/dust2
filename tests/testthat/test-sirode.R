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
  obj1 <- dust_system_create(sirode(), pars, n_particles = 1,
                            ode_control = ctl, deterministic = TRUE)
  dust_system_set_state_initial(obj1)
  dust_system_run_to_time(obj1, 9)
  d1 <- dust_system_internals(obj1)
  s1 <- dust_system_state(obj1)

  ## Here, I am get:
  ## observed: 4.1327704642489493
  ## expected: 4.353517

  ## so we're out by quite a bit, and I think that this is due to the
  ## interpolation being wrong, because I don't think we propagate
  ## that offset correctly!

  ## At the end of the step we have:
  ## c1 <- c(969.52426406888173, 20.026015297796405, 10.44972063332207,
  ##         20.475735931118471, 0.83294150528437161)
  ## c2 <- c(-7.4896420798189638, 3.612111020330893, 3.8775310594881027,
  ##         7.4896420798189922, 7.4896420798189949)
  ## c3 <- c(-0.60488124825999279, 0.27543100839365087, 0.32945023986631039, 
  ##         0.60488124825996437, 0.6048812482599617)
  ## c4 <- c(1.2097624965199856, -0.55086201678730173, -0.65890047973262078, 
  ##         -1.2097624965199287, -1.2097624965199234)
  ## c5 <- c(-0.00046681278128605895, 5.4702414467970872e-06,
  ##         0.00046134253984400403, 0.00046681278128605895,
  ##         0.00046681278128605895)

  ## It should be easy enough to create a deSolve version of this
  ## system that we can push around directly too.

  obj2 <- dust_system_create(sirode(), pars, n_particles = 1,
                            ode_control = ctl, deterministic = TRUE)
  dust_system_set_state_initial(obj2)
  invisible(dust_system_simulate(obj2, 0:9))
  d2 <- dust_system_internals(obj2)
  s2 <- dust_system_state(obj2)
})
