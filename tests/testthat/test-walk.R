test_that("can run simple walk model", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  expect_length(obj, 4)
  expect_type(obj[[1]], "externalptr")
  expect_equal(obj[[2]], 1)
  expect_false(obj[[3]])
  expect_null(obj[[4]])

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
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
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
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
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
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_set_state_initial(ptr))
  expect_equal(dust2_cpu_walk_state(ptr),
               rep(0, 10))
})


test_that("can run deterministically", {
  pars <- list(sd = 1)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, TRUE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_run_steps(ptr, 3))
  expect_equal(dust2_cpu_walk_state(ptr),
               rep(0, 10))
})


test_that("Allow fractional dt", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 0.5, 10, 0, 42, FALSE)
  ptr <- obj[[1]]
  expect_equal(dust2_cpu_walk_time(ptr), 0)
  expect_null(dust2_cpu_walk_run_steps(ptr, 3))
  expect_equal(dust2_cpu_walk_time(ptr), 1.5)
  expect_null(dust2_cpu_walk_run_steps(ptr, 5))
  expect_equal(dust2_cpu_walk_time(ptr), 4)
})


test_that("provided dt is reasonable", {
  pars <- list(sd = 1, random_initial = TRUE)
  expect_error(
    dust2_cpu_walk_alloc(pars, 0, 0, 10, 0, 42, FALSE),
    "Expected 'dt' to be greater than 0")
  expect_error(
    dust2_cpu_walk_alloc(pars, 0, -1, 10, 0, 42, FALSE),
    "Expected 'dt' to be greater than 0")
  expect_error(
    dust2_cpu_walk_alloc(pars, 0, 1.5, 10, 0, 42, FALSE),
    "Expected 'dt' to be at most 1")
  expect_error(
    dust2_cpu_walk_alloc(pars, 0, sqrt(2) / 2, 10, 0, 42, FALSE),
    "Expected 'dt' to be the inverse of an integer")
})


test_that("time starts as an integer", {
  pars <- list(sd = 1, random_initial = TRUE)
  expect_error(
    dust2_cpu_walk_alloc(pars, 1.5, 1, 10, 0, 42, FALSE),
    "Expected 'time' to be integer-like")
})


## Not a test at all of the model, but of the error handling in
## helpers.hpp
test_that("validate inputs", {
  pars <- list(sd = 1, random_initial = TRUE)
  expect_error(
    dust2_cpu_walk_alloc(pars, 0:3, 1, 10, 0, 42, FALSE),
    "'time' must be a scalar")

  expect_identical(
    dust2_cpu_walk_time(
      dust2_cpu_walk_alloc(pars, 5L, 1, 10, 0, 42, FALSE)[[1]]),
    5.0)
  expect_identical(
    dust2_cpu_walk_time(
      dust2_cpu_walk_alloc(pars, 5, 1, 10, 0, 42, FALSE)[[1]]),
    5.0)
  expect_error(
    dust2_cpu_walk_alloc(pars, "5", 1, 10, 0, 42, FALSE),
    "'time' must be scalar numeric")

  expect_identical(
    dust2_cpu_walk_state(
      dust2_cpu_walk_alloc(pars, 5, 1, 10, 0, 42, FALSE)[[1]]),
    rep(0, 10))
  expect_identical(
    dust2_cpu_walk_state(
      dust2_cpu_walk_alloc(pars, 5, 1, 10L, 0, 42, FALSE)[[1]]),
    rep(0, 10))
  expect_error(
    dust2_cpu_walk_alloc(pars, 5, 1, 9.5, 0, 42, FALSE),
    "'n_particles' must be integer-like")
  expect_error(
    dust2_cpu_walk_alloc(pars, 5, 1, "10", 0, 42, FALSE),
    "'n_particles' must be scalar integer")
  expect_error(
    dust2_cpu_walk_alloc(pars, 5, 1, -5, 0, 42, FALSE),
    "'n_particles' must be non-negative")

  expect_error(
    dust2_cpu_walk_alloc(pars, 5, 1, 10, 0, 42, 1),
    "'deterministic' must be scalar logical")
})


test_that("can initialise multiple groups with different parameter sets", {
  pars <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]
  expect_equal(dust2_cpu_walk_state(ptr), rep(0, 40))

  expect_null(dust2_cpu_walk_run_steps(ptr, 3))
  s <- dust2_cpu_walk_state(ptr)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 40)
  expect_equal(s, colSums(r$normal(3, 0, 1)) * rep(1:4, each = 10))
  expect_equal(dust2_cpu_walk_time(ptr), 3)
})


test_that("return names passed in with groups", {
  pars <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  names(pars) <- letters[1:4]
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 4, 42, FALSE)
  expect_true(obj[[3]])
  expect_equal(obj[[4]], letters[1:4])
})


test_that("can create multi-state walk model", {
  pars <- list(len = 3, sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  expect_equal(obj[[2]], 3)
  ptr <- obj[[1]]
  expect_equal(dust2_cpu_walk_state(ptr), rep(0, 30))
  expect_null(dust2_cpu_walk_set_state_initial(ptr))

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  s0 <- dust2_cpu_walk_state(ptr)
  expect_equal(s0, c(r$normal(3, 0, 1)))

  expect_null(dust2_cpu_walk_run_steps(ptr, 5))
  s1 <- dust2_cpu_walk_state(ptr)

  cmp <- r$normal(3 * 5, 0, 1)
  expect_equal(s1, s0 + c(apply(array(cmp, c(3, 5, 10)), c(1, 3), sum)))
})


test_that("prevent models from having different state lengths", {
  pars <- list(list(sd = 1, len = 2),
               list(sd = 2, len = 2),
               list(sd = 3, len = 3))
  expect_error(
    dust2_cpu_walk_alloc(pars, 0, 1, 10, 3, 42, FALSE),
    "Expected state length for group 3 to be 2, but it was 3")
})


test_that("require that parameter length matches requested number of groups", {
  pars <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  expect_error(
    dust2_cpu_walk_alloc(pars, 0, 1, 10, 3, 42, FALSE),
    "Expected 'pars' to have length 3 to match 'n_groups'")
})


test_that("can set time", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]
  expect_equal(dust2_cpu_walk_time(ptr), 0)
  expect_null(dust2_cpu_walk_set_time(ptr, 4))
  expect_equal(dust2_cpu_walk_time(ptr), 4)
  expect_null(dust2_cpu_walk_set_time(ptr, 0))
  expect_equal(dust2_cpu_walk_time(ptr), 0)
  expect_error(dust2_cpu_walk_set_time(ptr, 0.5),
               "Expected 'time' to be integer-like")
})


test_that("can update parameters", {
  pars1 <- list(sd = 1, random_initial = TRUE)
  pars2 <- list(sd = 10)
  obj <- dust2_cpu_walk_alloc(pars1, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]

  expect_null(dust2_cpu_walk_run_steps(ptr, 1))
  s1 <- dust2_cpu_walk_state(ptr)

  expect_null(dust2_cpu_walk_update_pars(ptr, pars2, FALSE))
  expect_null(dust2_cpu_walk_run_steps(ptr, 1))
  s2 <- dust2_cpu_walk_state(ptr)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  expect_equal(s1, drop(r$normal(1, 0, 1)))
  expect_equal(s2, s1 + drop(r$normal(1, 0, 10)))
})


test_that("can update parameters for grouped models", {
  pars1 <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  pars2 <- lapply(1:4, function(sd) list(sd = 10 * sd))

  obj <- dust2_cpu_walk_alloc(pars1, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]

  expect_null(dust2_cpu_walk_run_steps(ptr, 1))
  s1 <- dust2_cpu_walk_state(ptr)

  expect_null(dust2_cpu_walk_update_pars(ptr, pars2, TRUE))
  expect_null(dust2_cpu_walk_run_steps(ptr, 1))
  s2 <- dust2_cpu_walk_state(ptr)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 40)
  expect_equal(s1, drop(r$normal(1, 0, 1)) * rep(1:4, each = 10))
  expect_equal(s2, s1 + drop(r$normal(1, 0, 10)) * rep(1:4, each = 10))
})


test_that("can update parameters for grouped models", {
  pars1 <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  pars2 <- lapply(1:5, function(sd) list(sd = 10 * sd))
  obj <- dust2_cpu_walk_alloc(pars1, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]
  expect_error(dust2_cpu_walk_update_pars(ptr, pars2, TRUE),
               "Expected 'pars' to have length 4 to match 'n_groups'");
})
