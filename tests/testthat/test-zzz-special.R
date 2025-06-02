test_that("can run a system with special variables", {
  gen <- dust_compile("examples/special.cpp", quiet = TRUE, debug = TRUE)

  ## Initial conditions that make checking things easy; three
  ## particles starting from different points with qualitatively
  ## different initial values of the special rate variable.
  y0 <- rbind(c(0, 2, 4), c(-1, 0, 1))

  n_time <- 5
  n_particles <- 3

  sys <- dust_system_create(gen, list(), n_particles = n_particles, dt = 1)
  dust_system_set_state(sys, y0)
  t <- seq(0, n_time)
  y <- dust_system_simulate(sys, t)

  i_b <- dust_unpack_index(sys)$b

  sys0 <- dust_system_create(gen, list(), n_particles = n_particles)
  for (i in seq_len(n_time)) {
    dust_system_set_state(sys0, y[, , i])
    dust_system_run_to_time(sys0, t[[i + 1]])
    y0 <- dust_system_state(sys0)
    expect_equal(y0[-i_b, ], y[-i_b, , i + 1])
  }
})
