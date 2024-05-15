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

  expect_equal(dust2_cpu_walk_state(ptr, FALSE), matrix(0, 1, 10))
  expect_equal(dust2_cpu_walk_time(ptr), 0)

  expect_null(dust2_cpu_walk_run_steps(ptr, 3))
  s <- dust2_cpu_walk_state(ptr, FALSE)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  expect_equal(s, rbind(colSums(r$normal(3, 0, 1))))
  expect_equal(dust2_cpu_walk_time(ptr), 3)
})


test_that("can set model state from a vector", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]
  s <- rbind(runif(10))
  expect_null(dust2_cpu_walk_set_state(ptr, s, FALSE))
  expect_equal(dust2_cpu_walk_state(ptr, FALSE), s)

  expect_null(dust2_cpu_walk_run_steps(ptr, 3))

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  expect_equal(dust2_cpu_walk_state(ptr, FALSE),
               colSums(r$normal(3, 0, 1)) + s)
})


test_that("can set model state from initial conditions", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_set_state_initial(ptr))
  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  expect_equal(dust2_cpu_walk_state(ptr, FALSE),
               r$normal(1, 0, 1))
})


## Not 100% sure we'll keep this, but we will see.  This is only a
## feature of the walk model and not dust though.
test_that("can set model state from initial conditions with empty version", {
  pars <- list(sd = 1)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_set_state_initial(ptr))
  expect_equal(dust2_cpu_walk_state(ptr, FALSE),
               matrix(0, 1, 10))
})


test_that("can run deterministically", {
  pars <- list(sd = 1)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, TRUE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_run_steps(ptr, 3))
  expect_equal(dust2_cpu_walk_state(ptr, FALSE),
               matrix(0, 1, 10))
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
      dust2_cpu_walk_alloc(pars, 5, 1, 10, 0, 42, FALSE)[[1]], FALSE),
    matrix(0, 1, 10))
  expect_identical(
    dust2_cpu_walk_state(
      dust2_cpu_walk_alloc(pars, 5, 1, 10L, 0, 42, FALSE)[[1]], FALSE),
    matrix(0, 1, 10))
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
  expect_equal(dust2_cpu_walk_state(ptr, TRUE), array(0, c(1, 10, 4)))

  expect_null(dust2_cpu_walk_run_steps(ptr, 3))
  s <- dust2_cpu_walk_state(ptr, TRUE)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 40)
  expect_equal(
    s,
    array(colSums(r$normal(3, 0, 1)) * rep(1:4, each = 10), c(1, 10, 4)))
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
  expect_equal(dust2_cpu_walk_state(ptr, FALSE), matrix(0, 3, 10))
  expect_equal(dust2_cpu_walk_state(ptr, TRUE), array(0, c(3, 10, 1)))
  expect_null(dust2_cpu_walk_set_state_initial(ptr))

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  s0 <- dust2_cpu_walk_state(ptr, FALSE)
  expect_equal(s0, r$normal(3, 0, 1))

  expect_null(dust2_cpu_walk_run_steps(ptr, 5))
  s1 <- dust2_cpu_walk_state(ptr, FALSE)

  cmp <- r$normal(3 * 5, 0, 1)
  expect_equal(s1, s0 + apply(array(cmp, c(3, 5, 10)), c(1, 3), sum))
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
  s1 <- dust2_cpu_walk_state(ptr, FALSE)

  expect_null(dust2_cpu_walk_update_pars(ptr, pars2, FALSE))
  expect_null(dust2_cpu_walk_run_steps(ptr, 1))
  s2 <- dust2_cpu_walk_state(ptr, FALSE)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 10)
  expect_equal(s1, r$normal(1, 0, 1))
  expect_equal(s2, s1 + r$normal(1, 0, 10))
})


test_that("can update parameters for grouped models", {
  pars1 <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  pars2 <- lapply(1:4, function(sd) list(sd = 10 * sd))

  obj <- dust2_cpu_walk_alloc(pars1, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]

  expect_null(dust2_cpu_walk_run_steps(ptr, 1))
  s1 <- dust2_cpu_walk_state(ptr, FALSE)

  expect_null(dust2_cpu_walk_update_pars(ptr, pars2, TRUE))
  expect_null(dust2_cpu_walk_run_steps(ptr, 1))
  s2 <- dust2_cpu_walk_state(ptr, FALSE)

  r <- mcstate2::mcstate_rng$new(seed = 42, n_streams = 40)
  expect_equal(s1, r$normal(1, 0, 1) * rep(1:4, each = 10))
  expect_equal(s2, s1 + r$normal(1, 0, 10) * rep(1:4, each = 10))
})


test_that("params must be same length to update", {
  pars1 <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  pars2 <- lapply(1:5, function(sd) list(sd = 10 * sd))
  obj <- dust2_cpu_walk_alloc(pars1, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]
  expect_error(dust2_cpu_walk_update_pars(ptr, pars2, TRUE),
               "Expected 'pars' to have length 4 to match 'n_groups'")
})


test_that("can set state where n_state > 1", {
  pars <- list(len = 3, sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]

  ## One state per particle:
  s <- matrix(runif(30), 3, 10)
  expect_null(dust2_cpu_walk_set_state(ptr, s, FALSE))
  expect_equal(dust2_cpu_walk_state(ptr, FALSE), s)

  ## Shared state:
  s <- matrix(runif(3), 3, 1)
  expect_null(dust2_cpu_walk_set_state(ptr, s, FALSE))
  expect_equal(dust2_cpu_walk_state(ptr, FALSE), s[, rep(1, 10)])

  ## Appropriate errors:
  expect_error(
    dust2_cpu_walk_set_state(ptr, matrix(0, 2, 10), FALSE),
    "Expected the first dimension of 'state' to have size 3")
  expect_error(
    dust2_cpu_walk_set_state(ptr, matrix(0, 3, 15), FALSE),
    "Expected the second dimension of 'state' to have size 10 or 1")
  expect_error(
    dust2_cpu_walk_set_state(ptr, array(0, c(3, 10, 1)), FALSE),
    "Expected 'state' to be a 2d array")
})


test_that("can set state where n_state > 1 and groups are present", {
  pars <- lapply(1:4, function(i) list(len = 3, sd = i))
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]

  ## One state per particle:
  s <- array(runif(3 * 10 * 4), c(3, 10, 4))
  s <- array(as.numeric(seq_len(3 * 10 * 4)), c(3, 10, 4))
  expect_null(dust2_cpu_walk_set_state(ptr, s, TRUE))
  expect_equal(dust2_cpu_walk_state(ptr, TRUE), s)

  ## Shared state within groups
  s <- array(runif(3 * 1 * 4), c(3, 1, 4))
  s <- array(as.numeric(seq_len(3 * 1 * 4)), c(3, 1, 4))
  expect_null(dust2_cpu_walk_set_state(ptr, s, TRUE))
  expect_equal(dust2_cpu_walk_state(ptr, TRUE), s[, rep(1, 10), ])

  ## Shared state across groups
  s <- array(runif(3 * 10 * 1), c(3, 10, 1))
  s <- array(as.numeric(seq_len(3 * 10 * 1)), c(3, 10, 1))
  expect_null(dust2_cpu_walk_set_state(ptr, s, TRUE))
  expect_equal(dust2_cpu_walk_state(ptr, TRUE), s[, , rep(1, 4)])

  ## Shared state across everything
  s <- array(runif(3 * 1 * 1), c(3, 1, 1))
  s <- array(as.numeric(seq_len(3 * 1 * 1)), c(3, 1, 1))
  expect_null(dust2_cpu_walk_set_state(ptr, s, TRUE))
  expect_equal(dust2_cpu_walk_state(ptr, TRUE), s[, rep(1, 10), rep(1, 4)])

  ## Appropriate errors:
  expect_error(
    dust2_cpu_walk_set_state(ptr, array(0, c(2, 10, 4)), TRUE),
    "Expected the first dimension of 'state' to have size 3")
  expect_error(
    dust2_cpu_walk_set_state(ptr, array(0, c(3, 15, 4)), TRUE),
    "Expected the second dimension of 'state' to have size 10 or 1")
  expect_error(
    dust2_cpu_walk_set_state(ptr, array(0, c(3, 10, 2)), TRUE),
    "Expected the third dimension of 'state' to have size 4 or 1")
  expect_error(
    dust2_cpu_walk_set_state(ptr, array(0, c(3, 10)), TRUE),
    "Expected 'state' to be a 3d array")
})


## This one is for consistency with some of the bits above and allows
## setting state while ignoring the grouped structure.  This is most
## likely to be used when having one particle per group but it can be
## used in any context really.
test_that("can set ungrouped state into grouped model", {
  pars <- lapply(1:4, function(i) list(len = 3, sd = i))
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]

  ## One state per particle:
  s <- array(runif(30), c(3, 10 * 4))
  expect_null(dust2_cpu_walk_set_state(ptr, s, FALSE))
  expect_equal(dust2_cpu_walk_state(ptr, FALSE), s)
  expect_equal(dust2_cpu_walk_state(ptr, TRUE),
               array(s, c(3, 10, 4)))

  s <- array(runif(3), c(3, 1))
  expect_null(dust2_cpu_walk_set_state(ptr, s, FALSE))
  expect_equal(dust2_cpu_walk_state(ptr, FALSE), array(s, c(3, 40)))
  expect_equal(dust2_cpu_walk_state(ptr, TRUE),
               array(s, c(3, 10, 4)))
})


test_that("can reorder state", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_set_state_initial(ptr))
  s1 <- dust2_cpu_walk_state(ptr, FALSE)
  i <- sample(10, replace = TRUE)
  expect_null(dust2_cpu_walk_reorder(ptr, i))
  s2 <- dust2_cpu_walk_state(ptr, FALSE)
  expect_equal(s2, s1[, i, drop = FALSE])
})


test_that("can reorder state in a multiparameter model", {
  pars <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 4, 42, FALSE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_set_state_initial(ptr))
  s1 <- dust2_cpu_walk_state(ptr, TRUE)
  i <- replicate(4, sample(10, replace = TRUE))
  expect_null(dust2_cpu_walk_reorder(ptr, i))
  s2 <- dust2_cpu_walk_state(ptr, TRUE)
  expect_equal(
    s2,
    vapply(1:4, function(j) s1[, i[, j], j, drop = FALSE], matrix(0, 1, 10)))
})


test_that("validate inputs for reordering", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 0, 42, FALSE)
  ptr <- obj[[1]]
  expect_null(dust2_cpu_walk_set_state_initial(ptr))
  s1 <- dust2_cpu_walk_state(ptr, FALSE)
  expect_error(
    dust2_cpu_walk_reorder(ptr, c(1L, 2L, 3L)),
    "Expected an index of length 10")
  expect_error(
    dust2_cpu_walk_reorder(ptr, seq_len(10) + 1L),
    "Expected 'index' values to lie in [1, 10]",
    fixed = TRUE)
})


test_that("can run walk model to time", {
  pars <- list(sd = 1)
  obj1 <- dust2_cpu_walk_alloc(pars, 0, 0.25, 10, 0, 42, FALSE)
  obj2 <- dust2_cpu_walk_alloc(pars, 0, 0.25, 10, 0, 42, FALSE)
  ptr1 <- obj1[[1]]
  ptr2 <- obj2[[1]]

  dust2_cpu_walk_set_state_initial(ptr1)
  expect_null(dust2_cpu_walk_run_steps(ptr1, 40))
  expect_equal(dust2_cpu_walk_time(ptr1), 10)
  s1 <- dust2_cpu_walk_state(ptr1, FALSE)

  dust2_cpu_walk_set_state_initial(ptr2)
  expect_null(dust2_cpu_walk_run_to_time(ptr2, 10))
  expect_equal(dust2_cpu_walk_time(ptr2), 10)
  expect_equal(dust2_cpu_walk_state(ptr2, FALSE), s1)
})
