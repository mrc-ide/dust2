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
