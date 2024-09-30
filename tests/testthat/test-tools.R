test_that("can unpack state from systems with a single particle", {
  sys <- dust_system_create(sir(), list(), n_particles = 1)
  dust_system_set_state_initial(sys)
  s <- dust_system_state(sys)
  s2 <- dust_unpack_state(sys, s)
  expect_equal(s2, sys$packer_state$unpack(s))
  expect_equal(s2, list(S = 990, I = 10, R = 0, cases_cumul = 0, cases_inc = 0))
})


test_that("can unpack state from systems with several particles", {
  sys <- dust_system_create(sir(), list(), n_particles = 10)
  dust_system_set_state_initial(sys)
  s <- dust_system_state(sys)
  s2 <- dust_unpack_state(sys, s)
  expect_equal(s2, sys$packer_state$unpack(s))
  expect_equal(lengths(s2, FALSE), rep(10, 5))
})


test_that("can't get unpacker from filter before running", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_filter_create(sir(), time_start, data, n_particles = 10)
  expect_error(
    dust_unpack_state(obj, numeric()),
    "Packer is not yet ready")
})


test_that("can't get unpacker from unknown object", {
  expect_error(
    dust_unpack_state(NULL, numeric()),
    "Expected 'obj' to be a 'dust_system' or a 'dust_likelihood'")
})
