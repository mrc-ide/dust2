test_that("can run simple walk model", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 42, FALSE)
  expect_type(obj[[1]], "externalptr")
  expect_equal(obj[[2]], 1)

  ptr <- obj[[1]]
  expect_type(dust2_cpu_walk_rng_state(ptr), "raw")
  expect_length(dust2_cpu_walk_rng_state(ptr), 32 * 10)

  expect_equal(dust2_cpu_walk_state(ptr), rep(0, 10))
  expect_equal(dust2_cpu_walk_time(ptr), 0)

  expect_null(dust2_cpu_walk_run_steps(ptr, 3))
  s <- dust2_cpu_walk_state(ptr)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  expect_equal(s, colSums(r$normal(3, 0, 1)))
  expect_equal(dust2_cpu_walk_time(ptr), 3)
})
