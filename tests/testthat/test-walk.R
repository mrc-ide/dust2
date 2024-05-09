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


test_that("can set model state from a vector", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 42, FALSE)
  ptr <- obj[[1]]
  s <- runif(10)
  expect_null(dust2_cpu_walk_set_state(ptr, s))
  expect_equal(dust2_cpu_walk_state(ptr), s)

  expect_null(dust2_cpu_walk_run_steps(ptr, 3))

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  expect_equal(dust2_cpu_walk_state(ptr),
               colSums(r$normal(3, 0, 1)) + s)
})


test_that("can set model state from initial conditions", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 42, FALSE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_set_state_initial(ptr))
  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  expect_equal(dust2_cpu_walk_state(ptr),
               drop(r$normal(1, 0, 1)))
})


## Not 100% sure we'll keep this, but we will see.  This is only a
## feature of the walk model and not dust though.
test_that("can set model state from initial conditions with empty version", {
  pars <- list(sd = 1)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 42, FALSE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_set_state_initial(ptr))
  expect_equal(dust2_cpu_walk_state(ptr),
               rep(0, 10))
})


test_that("can run deterministically", {
  pars <- list(sd = 1)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 42, TRUE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_run_steps(ptr, 3))
  expect_equal(dust2_cpu_walk_state(ptr),
               rep(0, 10))
})


test_that("require that dt is 1 for now", {
  pars <- list(sd = 1, random_initial = TRUE)
  expect_error(
    dust2_cpu_walk_alloc(pars, 0, 0.5, 10, 42, FALSE),
    "Requiring dt = 1 for now",
    fixed = TRUE)
})


## Not a test at all of the model, but of the error handling in
## helpers.hpp
test_that("validate inputs", {
  pars <- list(sd = 1, random_initial = TRUE)
  expect_error(
    dust2_cpu_walk_alloc(pars, 0:3, 1, 10, 42, FALSE),
    "'time' must be a scalar")

  expect_identical(
    dust2_cpu_walk_time(dust2_cpu_walk_alloc(pars, 5L, 1, 10, 42, FALSE)[[1]]),
    5.0)
  expect_identical(
    dust2_cpu_walk_time(dust2_cpu_walk_alloc(pars, 5, 1, 10, 42, FALSE)[[1]]),
    5.0)
  expect_error(
    dust2_cpu_walk_alloc(pars, "5", 1, 10, 42, FALSE),
    "'time' must be scalar numeric")

  expect_identical(
    dust2_cpu_walk_state(dust2_cpu_walk_alloc(pars, 5, 1, 10, 42, FALSE)[[1]]),
    rep(0, 10))
  expect_identical(
    dust2_cpu_walk_state(dust2_cpu_walk_alloc(pars, 5, 1, 10L, 42, FALSE)[[1]]),
    rep(0, 10))
  expect_error(
    dust2_cpu_walk_alloc(pars, 5, 1, 9.5, 42, FALSE),
    "'n_particles' must be integer-like")
  expect_error(
    dust2_cpu_walk_alloc(pars, 5, 1, "10", 42, FALSE),
    "'n_particles' must be scalar integer")
  expect_error(
    dust2_cpu_walk_alloc(pars, 5, 1, -5, 42, FALSE),
    "'n_particles' must be non-negative")

  expect_error(
    dust2_cpu_walk_alloc(pars, 5, 1, 10, 42, 1),
    "'deterministic' must be scalar logical")
})
