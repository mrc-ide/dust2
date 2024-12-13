test_that("can validate that 'dt' is reasonable", {
  expect_no_error(check_dt(1))
  expect_no_error(check_dt(1 / 5))
  expect_no_error(check_dt(2))
  expect_no_error(check_dt(NULL, allow_null = TRUE))
  expect_error(check_dt(-1, name = "dt"),
               "Expected 'dt' to be greater than 0")
  expect_error(check_dt(2.5, name = "dt"),
               "Expected 'dt' to be an integer, if greater than 1")
  expect_error(check_dt(1 / 3.5, name = "dt"),
               "Expected 'dt' to be the inverse of an integer")
  expect_error(check_dt(NULL, name = "dt"),
               "'dt' must be a scalar")
  expect_error(check_dt(Inf, name = "dt"),
               "Expected 'dt' to be an integer, if greater than 1")
})


test_that("check index", {
  expect_no_error(check_index(NULL))
  expect_no_error(check_index(1:4))
  idx <- c(1, 2, 3.4)
  expect_error(check_index(idx),
               "Expected 'idx' to be integer")
  idx <- c(1, 2, -3)
  expect_error(check_index(idx),
               "All elements of 'idx' must be at least 1")
})

test_that("can apply stochastic updates by setting dt", {
  sys <- dust_system_create(malaria(), list(), n_particles = 10, dt = 10)
  dust_system_set_state_initial(sys)
  t <- seq(0, 20)
  y <- dust_system_simulate(sys, t)
  i_beta <- dust_unpack_index(sys)$beta
  beta <- y[i_beta, , ]
  
  ## Equivalence to a deterinistic model:
  sys0 <- dust_system_create(malaria(), list(), n_particles = 10)
  for (i in seq_len(length(t) - 1)) {
    dust_system_set_state(sys0, y[, , i])
    dust_system_run_to_time(sys0, t[[i + 1]])
    y0 <- dust_system_state(sys0)
    #expect_equal(y0[-i_beta, ], y[-i_beta, , i + 1])
  }
})