test_that("can run model without any stochastic updates", {
  sys <- dust_system_create(malaria(), list(), n_particles = 10)
  dust_system_set_state_initial(sys)
  t <- seq(0, 20)
  y <- dust_system_simulate(sys, t)
  yy <- dust_unpack_state(sys, y)
  expect_equal(y, y[, rep(1, 10), , drop = FALSE])
  expect_equal(drop(yy$beta), array(log(10 / 9), c(10, 21)))
})


test_that("can apply stochastic updates by setting dt", {
  sys <- dust_system_create(malaria(), list(), n_particles = 10, dt = 1)
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
    expect_equal(y0[-i_beta, ], y[-i_beta, , i + 1])
  }
})


test_that("can print a mixed generator", {
  res <- evaluate_promise(withVisible(print(malaria)))
  expect_mapequal(res$result, list(value = malaria, visible = FALSE))
  expect_match(res$messages, "<dust_system_generator: malaria>",
               fixed = TRUE, all = FALSE)
  expect_match(
    res$messages,
    "This system runs in both continuous time and discrete time",
    all = FALSE)
})


test_that("Can compare to data", {
  sys <- dust_system_create(malaria(), list(), n_particles = 10, dt = 1)
  dust_system_set_state_initial(sys)
  dust_system_run_to_time(sys, 10)
  s <- dust_unpack_state(sys, dust_system_state(sys))
  d <- list(tested = 4, positive = 2)
  expect_equal(dust_system_compare_data(sys, d),
               dbinom(d$positive, d$tested, s$Ih, log = TRUE))
})
