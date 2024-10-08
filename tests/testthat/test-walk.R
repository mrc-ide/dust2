test_that("can run simple walk system", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_s3_class(obj, "dust_system")

  expect_type(dust_system_rng_state(obj), "raw")
  expect_length(dust_system_rng_state(obj), 32 * 10)

  expect_equal(dust_system_state(obj), matrix(0, 1, 10))
  expect_equal(dust_system_time(obj), 0)

  expect_null(dust_system_run_to_time(obj, 3))
  s <- dust_system_state(obj)

  r <- monty::monty_rng$new(seed = 42, n_streams = 10)
  expect_equal(s, rbind(colSums(r$normal(3, 0, 1))))
  expect_equal(dust_system_time(obj), 3)

  expect_identical(dim(obj), c(1L, 10L))
})


test_that("can set system state from a vector", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  s <- rbind(runif(10))
  expect_null(dust_system_set_state(obj, s))
  expect_equal(dust_system_state(obj), s)

  expect_null(dust_system_run_to_time(obj, 3))

  r <- monty::monty_rng$new(seed = 42, n_streams = 10)
  expect_equal(dust_system_state(obj),
               colSums(r$normal(3, 0, 1)) + s)
})


test_that("can set system state from initial conditions", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_null(dust_system_set_state_initial(obj))
  r <- monty::monty_rng$new(seed = 42, n_streams = 10)
  expect_equal(dust_system_state(obj),
               r$normal(1, 0, 1))
})


## Not 100% sure we'll keep this, but we will see.  This is only a
## feature of the walk system and not dust though.
test_that("can set system state from initial conditions with empty version", {
  pars <- list(sd = 1)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_null(dust_system_set_state_initial(obj))
  expect_equal(dust_system_state(obj),
               matrix(0, 1, 10))
})


test_that("can run deterministically", {
  pars <- list(sd = 1)
  obj <- dust_system_create(walk(), pars, n_particles = 10,
                            deterministic = TRUE)
  expect_null(dust_system_run_to_time(obj, 3))
  expect_equal(dust_system_state(obj),
               matrix(0, 1, 10))
})


test_that("Allow running with dt", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, dt = 0.5, n_particles = 10, seed = 42)
  expect_equal(dust_system_time(obj), 0)
  expect_null(dust_system_run_to_time(obj, 2))
  expect_equal(dust_system_time(obj), 2)
})


test_that("provided dt is reasonable", {
  pars <- list(sd = 1, random_initial = TRUE)
  expect_error(
    dust_system_create(walk(), pars, n_particles = 10, dt = 0),
    "Expected 'dt' to be greater than 0")
  expect_error(
    dust_system_create(walk(), pars, n_particles = 10, dt = -1),
    "Expected 'dt' to be greater than 0")
  expect_error(
    dust_system_create(walk(), pars, n_particles = 10, dt = 1.5),
    "Expected 'dt' to be at most 1")
  expect_error(
    dust_system_create(walk(), pars, n_particles = 10, dt = sqrt(2) / 2),
    "Expected 'dt' to be the inverse of an integer")
})


test_that("time starts as an integer when dt is 1", {
  pars <- list(sd = 1, random_initial = TRUE)
  expect_error(
    dust_system_create(walk(), pars, n_particles = 10, time = 1.5),
    "'time' must be integer-like, because 'dt' is 1")
})


test_that("time aligns to grid when dt is less than 1", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, dt = 0.1,
                            time = 1.5)
  expect_equal(dust_system_time(obj), 1.5)
  expect_error(
    dust_system_create(walk(), pars, n_particles = 10, dt = 0.25, time = 1.1),
    "'time' must be a multiple of 'dt' (0.25)",
    fixed = TRUE)
})


## Not a test at all of the system, but of the error handling in
## helpers.hpp
test_that("validate inputs", {
  pars <- list(sd = 1, random_initial = TRUE)
  expect_error(
    dust_system_create(walk(), pars, time = 1:3, n_particles = 10),
    "'time' must be a scalar")

  expect_identical(
    dust_system_time(
      dust_system_create(walk(), pars, time = 5L, n_particles = 10)),
    5.0)
  expect_identical(
    dust_system_time(
      dust_system_create(walk(), pars, time = 5, n_particles = 10)),
    5.0)
  expect_error(
    dust_system_create(walk(), pars, time = "5", n_particles = 10),
    "Expected 'time' to be numeric")

  expect_identical(
    dust_system_state(
      dust_system_create(walk(), pars, n_particles = 10)),
    matrix(0, 1, 10))
  expect_identical(
    dust_system_state(
      dust_system_create(walk(), pars, n_particles = 10L)),
    matrix(0, 1, 10))
  expect_error(
    dust_system_create(walk(), pars, n_particles = 9.5),
    "Expected 'n_particles' to be integer")
  expect_error(
    dust_system_create(walk(), pars, n_particles = "10"),
    "Expected 'n_particles' to be integer")
  expect_error(
    dust_system_create(walk(), pars, n_particles = -5),
    "'n_particles' must be at least 1")

  expect_error(
    dust_system_create(walk(), pars, n_particles = 10, deterministic = 1),
    "'deterministic' must be scalar logical")

  expect_error(
    dust_system_create(walk(), list(), n_particles = 10),
    "A value is expected for 'sd'")
})


test_that("can initialise multiple groups with different parameter sets", {
  pars <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  obj <- dust_system_create(walk(), pars, n_particles = 10, n_groups = 4,
                           seed = 42)
  expect_equal(dust_system_state(obj), array(0, c(1, 10, 4)))

  expect_null(dust_system_run_to_time(obj, 3))
  s <- dust_system_state(obj)

  r <- monty::monty_rng$new(seed = 42, n_streams = 40)
  expect_equal(
    s,
    array(colSums(r$normal(3, 0, 1)) * rep(1:4, each = 10), c(1, 10, 4)))
  expect_equal(dust_system_time(obj), 3)

  expect_identical(dim(obj), c(1L, 10L, 4L))
})


test_that("can create multi-state walk system", {
  pars <- list(len = 3, sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_equal(obj$n_state, 3)
  expect_equal(dust_system_state(obj), array(0, c(3, 10)))
  expect_null(dust_system_set_state_initial(obj))

  r <- monty::monty_rng$new(seed = 42, n_streams = 10)
  s0 <- dust_system_state(obj)
  expect_equal(s0, r$normal(3, 0, 1))

  expect_null(dust_system_run_to_time(obj, 5))
  s1 <- dust_system_state(obj)

  cmp <- r$normal(3 * 5, 0, 1)
  expect_equal(s1, s0 + apply(array(cmp, c(3, 5, 10)), c(1, 3), sum))
  expect_identical(dim(obj), c(3L, 10L))
})


test_that("prevent systems from having different state lengths", {
  pars <- list(list(sd = 1, len = 2),
               list(sd = 2, len = 2),
               list(sd = 3, len = 3))
  expect_error(
    dust_system_create(walk(), pars, n_particles = 10, n_groups = 3),
    paste("State length for group 3 was different to previous groups; total",
          "length was expected to be 2 but it was 3"))
})


test_that("require that parameter length matches requested number of groups", {
  pars <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  expect_error(
    dust_system_create(walk(), pars, n_particles = 10, n_groups = 3),
    "Expected 'pars' to have length 3 to match 'n_groups'")
})


test_that("can set time", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10)
  expect_equal(dust_system_time(obj), 0)
  expect_null(dust_system_set_time(obj, 4))
  expect_equal(dust_system_time(obj), 4)
  expect_null(dust_system_set_time(obj, 0))
  expect_equal(dust_system_time(obj), 0)
  expect_error(dust_system_set_time(obj, 0.5),
               "'time' must be integer-like, because 'dt' is 1")
})


test_that("can set fractional time", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, dt = 0.25)
  expect_error(
    dust_system_set_time(obj, 0.1),
    "'time' must be a multiple of 'dt' (0.25)",
    fixed = TRUE)
  dust_system_set_time(obj, 1.25)
  expect_equal(dust_system_time(obj), 1.25)
})


test_that("can update parameters", {
  pars1 <- list(sd = 1, random_initial = TRUE)
  pars2 <- list(sd = 10)
  obj <- dust_system_create(walk(), pars1, n_particles = 10, seed = 42)

  expect_null(dust_system_run_to_time(obj, 1))
  s1 <- dust_system_state(obj)

  expect_null(dust_system_update_pars(obj, pars2))
  expect_null(dust_system_run_to_time(obj, 2))
  s2 <- dust_system_state(obj)

  r <- monty::monty_rng$new(seed = 42, n_streams = 10)
  expect_equal(s1, r$normal(1, 0, 1))
  expect_equal(s2, s1 + r$normal(1, 0, 10))
})


test_that("can update parameters for grouped systems", {
  pars1 <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  pars2 <- lapply(1:4, function(sd) list(sd = 10 * sd))

  obj <- dust_system_create(walk(), pars1, n_particles = 10, n_groups = 4,
                           seed = 42)

  expect_null(dust_system_run_to_time(obj, 1))
  s1 <- dust_system_state(obj)

  expect_null(dust_system_update_pars(obj, pars2))
  expect_null(dust_system_run_to_time(obj, 2))
  s2 <- dust_system_state(obj)

  r <- monty::monty_rng$new(seed = 42, n_streams = 40)
  expect_equal(s1, array(r$normal(1, 0, 1) * rep(1:4, each = 10),
                         c(1, 10, 4)))
  expect_equal(s2, array(s1 + drop(r$normal(1, 0, 10)) * rep(1:4, each = 10),
                         c(1, 10, 4)))
})


test_that("params must be same length to update", {
  pars1 <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  pars2 <- lapply(1:5, function(sd) list(sd = 10 * sd))
  obj <- dust_system_create(walk(), pars1, n_particles = 10, n_groups = 4)
  expect_error(dust_system_update_pars(obj, pars2),
               "Expected 'pars' to have length 4 to match 'n_groups'")
})


test_that("can set state where n_state > 1", {
  pars <- list(len = 3, sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)

  ## One state per particle:
  s <- matrix(runif(30), 3, 10)
  expect_null(dust_system_set_state(obj, s))
  expect_equal(dust_system_state(obj), s)

  ## Shared state:
  s <- matrix(runif(3), 3, 1)
  expect_null(dust_system_set_state(obj, s))
  expect_equal(dust_system_state(obj), s[, rep(1, 10)])

  ## Appropriate errors:
  expect_error(
    dust_system_set_state(obj, matrix(0, 2, 10)),
    "Expected the first dimension of 'state' to have size 3")
  expect_error(
    dust_system_set_state(obj, matrix(0, 3, 15)),
    "Expected the second dimension of 'state' to have size 10 or 1")
  expect_error(
    dust_system_set_state(obj, array(0, c(3, 10, 1))),
    "Expected 'state' to be a 2d array")
})


test_that("can set state where n_state > 1 and groups are present", {
  pars <- lapply(1:4, function(i) list(len = 3, sd = i))
  obj <- dust_system_create(walk(), pars, n_particles = 10, n_groups = 4,
                           seed = 42)

  ## One state per particle:
  s <- array(runif(3 * 10 * 4), c(3, 10, 4))
  s <- array(as.numeric(seq_len(3 * 10 * 4)), c(3, 10, 4))
  expect_null(dust_system_set_state(obj, s))
  expect_equal(dust_system_state(obj), s)

  ## Shared state within groups
  s <- array(runif(3 * 1 * 4), c(3, 1, 4))
  s <- array(as.numeric(seq_len(3 * 1 * 4)), c(3, 1, 4))
  expect_null(dust_system_set_state(obj, s))
  expect_equal(dust_system_state(obj), s[, rep(1, 10), ])

  ## Shared state across groups
  s <- array(runif(3 * 10 * 1), c(3, 10, 1))
  s <- array(as.numeric(seq_len(3 * 10 * 1)), c(3, 10, 1))
  expect_null(dust_system_set_state(obj, s))
  expect_equal(dust_system_state(obj), s[, , rep(1, 4)])

  ## Shared state across everything
  s <- array(runif(3 * 1 * 1), c(3, 1, 1))
  s <- array(as.numeric(seq_len(3 * 1 * 1)), c(3, 1, 1))
  expect_null(dust_system_set_state(obj, s))
  expect_equal(dust_system_state(obj), s[, rep(1, 10), rep(1, 4)])

  ## Appropriate errors:
  expect_error(
    dust_system_set_state(obj, array(0, c(2, 10, 4))),
    "Expected the first dimension of 'state' to have size 3")
  expect_error(
    dust_system_set_state(obj, array(0, c(3, 15, 4))),
    "Expected the second dimension of 'state' to have size 10 or 1")
  expect_error(
    dust_system_set_state(obj, array(0, c(3, 10, 2))),
    "Expected the third dimension of 'state' to have size 4 or 1")
  expect_error(
    dust_system_set_state(obj, array(0, c(3, 10))),
    "Expected 'state' to be a 3d array")
})


test_that("can reorder state", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_null(dust_system_set_state_initial(obj))
  s1 <- dust_system_state(obj)
  i <- sample(10, replace = TRUE)
  expect_null(dust_system_reorder(obj, i))
  s2 <- dust_system_state(obj)
  expect_equal(s2, s1[, i, drop = FALSE])
})


test_that("can reorder state in a multiparameter system", {
  pars <- lapply(1:4, function(sd) list(sd = sd, random_initial = TRUE))
  obj <- dust_system_create(walk(), pars, n_particles = 10, n_groups = 4,
                           seed = 42)
  expect_null(dust_system_set_state_initial(obj))
  s1 <- dust_system_state(obj)
  i <- replicate(4, sample(10, replace = TRUE))
  expect_null(dust_system_reorder(obj, i))
  s2 <- dust_system_state(obj)
  expect_equal(
    s2,
    vapply(1:4, function(j) s1[, i[, j], j, drop = FALSE], matrix(0, 1, 10)))
})


test_that("validate inputs for reordering", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_null(dust_system_set_state_initial(obj))
  s1 <- dust_system_state(obj)
  expect_error(
    dust_system_reorder(obj, c(1L, 2L, 3L)),
    "Expected an index of length 10")
  expect_error(
    dust_system_reorder(obj, seq_len(10) + 1L),
    "Expected 'index' values to lie in [1, 10]",
    fixed = TRUE)
})


test_that("time must not be in the past", {
  pars <- list(sd = 1)
  obj <- dust_system_create(walk(), pars, dt = 0.25, n_particles = 10)

  dust_system_set_state_initial(obj)
  expect_error(
    dust_system_run_to_time(obj, -5),
    "Can't run to time -5.*, system already at time 0.*")
})


test_that("can simulate walk system", {
  pars <- lapply(1:4, function(sd) list(len = 3, sd = sd))
  obj1 <- dust_system_create(walk(), pars, n_particles = 10, n_groups = 4,
                            seed = 42)
  obj2 <- dust_system_create(walk(), pars, n_particles = 10, n_groups = 4,
                            seed = 42)

  res <- dust_system_simulate(obj1, 0:20)
  expect_equal(dim(res), c(3, 10, 4, 21))

  expect_equal(res[, , , 21], dust_system_state(obj1))
  expect_equal(dust_system_time(obj1), 20)

  expect_equal(res[, , , 1], dust_system_state(obj2))
  dust_system_run_to_time(obj2, 10)
  expect_equal(res[, , , 11], dust_system_state(obj2))
})


test_that("can simulate walk system with index", {
  pars <- lapply(1:4, function(sd) list(len = 6, sd = sd))
  obj1 <- dust_system_create(walk(), pars, n_particles = 10, n_groups = 4,
                            seed = 42)
  obj2 <- dust_system_create(walk(), pars, n_particles = 10, n_groups = 4,
                            seed = 42)

  res1 <- dust_system_simulate(obj1, 0:20, c(2, 4))
  res2 <- dust_system_simulate(obj2, 0:20)
  expect_equal(res1, res2[c(2, 4), , , ])
})


test_that("can validate times", {
  pars <- list(sd = 1)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_error(
    dust_system_simulate(obj, c(1, 2, 3.5)),
    "Values in 'times' must be integer-like, because 'dt' is 1",
    fixed = TRUE)
  expect_error(
    dust_system_simulate(obj, c(1, 2, 1)),
    "Values in 'times' must be increasing",
    fixed = TRUE)
  expect_error(
    dust_system_simulate(obj, NULL),
    "Expected 'times' to be numeric")
  expect_error(
    dust_system_simulate(obj, numeric()),
    "Expected at least one value in 'times'")
})


test_that("allow non-integer times but validate times align", {
  pars <- list(sd = 1)
  obj1 <- dust_system_create(walk(), pars, n_particles = 10, seed = 42,
                             dt = 0.25)
  obj2 <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_error(
    dust_system_simulate(obj1, c(1, 2, 3.1)),
    "Values in 'times' must be multiples of 'dt'")

  t <- seq(0, 1, by = 0.25)
  y1 <- dust_system_simulate(obj1, t)
  y2 <- dust_system_simulate(obj2, t * 4) / 4
  expect_equal(y1, y2)
})


test_that("can validate index values", {
  pars <- list(len = 3, sd = 1)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_error(
    dust_system_simulate(obj, 0:20, c("2", "4")),
    "Expected an integer vector for 'index_state'")
  expect_error(
    dust_system_simulate(obj, 0:20, c(2, 2.9)),
    paste("All values of 'index_state' must be integer-like, but",
          "'index_state[2]' was not"),
    fixed = TRUE)
  expect_error(
    dust_system_simulate(obj, 0:20, c(2, 20)),
    paste("All values of 'index_state' must be in [1, 3], but",
          "'index_state[2]' was 20"),
          fixed = TRUE)
  expect_error(
    dust_system_simulate(obj, 0:20, c(2, 3, -4)),
    paste("All values of 'index_state' must be in [1, 3], but",
          "'index_state[3]' was -4"),
    fixed = TRUE)
  expect_error(
    dust_system_simulate(obj, 0:20, integer()),
    "'index_state' must have nonzero length")
})


test_that("can't compare walk to data", {
  pars <- list(sd = 1)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_error(
    dust_system_compare_data(obj, list()),
    paste("Can't compare against data; the 'walk' system does not have",
          "'compare_data'"))
})


test_that("can't create filter with walk system", {
  pars <- list(sd = 1)
  time_start <- 0
  data <- list(time = 1:4, d = 1:4)
  expect_error(
    dust_unfilter_create(walk(), pars, time_start, data),
    "Can't create unfilter; the 'walk' system does not have 'compare_data'")
  expect_error(
    dust_filter_create(walk(), pars, time_start, data, 10),
    "Can't create filter; the 'walk' system does not have 'compare_data'")
})


test_that("can set rng state", {
  pars <- list(sd = 1, random_initial = TRUE)

  obj1 <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  obj2 <- dust_system_create(walk(), pars, n_particles = 10, seed = 43)

  expect_false(identical(dust_system_rng_state(obj1),
                         dust_system_rng_state(obj2)))

  expect_null(dust_system_set_rng_state(obj2, dust_system_rng_state(obj1)))
  expect_identical(dust_system_rng_state(obj1), dust_system_rng_state(obj2))
  dust_system_run_to_time(obj1, 10)
  dust_system_run_to_time(obj2, 10)
  expect_identical(dust_system_state(obj1), dust_system_state(obj2))
})


test_that("can validate rng state on setting", {
  pars <- list(sd = 1)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  s <- dust_system_rng_state(obj)
  expect_error(dust_system_set_rng_state(obj, NULL),
               "Expected a raw vector for 'rng_state'")
  expect_error(dust_system_set_rng_state(obj, s[-1]),
               "Incorrect length for 'rng_state'; expected 320 but given 319")
})


test_that("reject ode_control as an argument when creating walk model", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_error(
    dust_system_create(walk(), pars, ode_control = dust_ode_control(),
                       n_particles = 1),
    "Can't use 'ode_control' with discrete-time systems")
})


test_that("discrete time models have no internals", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 10, seed = 42)
  expect_null(dust_system_run_to_time(obj, 3))
  expect_null(dust_system_internals(obj))
})


test_that("can extract subset of state from an ungrouped system", {
  set.seed(1)
  pars <- list(len = 5, sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 7)
  dust_system_set_state_initial(obj)
  m <- dust_system_state(obj)
  expect_equal(dim(m), c(5, 7))

  expect_equal(
    dust_system_state(obj, index_particle = c(2, 5, 6)),
    m[, c(2, 5, 6)])
  expect_equal(
    dust_system_state(obj, index_state = c(1, 2, 4)),
    m[c(1, 2, 4), ])
  expect_equal(
    dust_system_state(obj, index_particle = c(4, 7), index_state = c(1, 2, 4)),
    m[c(1, 2, 4), c(4, 7)])

  expect_error(
    dust_system_state(obj, index_particle = integer(0)),
    "'index_particle' must have nonzero length")
  expect_error(
    dust_system_state(obj, index_state = integer(0)),
    "'index_state' must have nonzero length")

  expect_error(
    dust_system_state(obj, index_group = 1),
    "Can't provide 'index_group' for a non-grouped system")
  expect_error(
    dust_system_state(obj, index_group = 1:2),
    "Can't provide 'index_group' for a non-grouped system")
})


test_that("can extract subset of state from a grouped system", {
  set.seed(1)
  pars <- rep(list(list(len = 5, sd = 1, random_initial = TRUE)), 3)
  obj <- dust_system_create(walk(), pars, n_particles = 7, n_groups = 3)
  dust_system_set_state_initial(obj)
  m <- dust_system_state(obj)
  expect_equal(dim(m), c(5, 7, 3))

  expect_equal(
    dust_system_state(obj, index_particle = c(2, 5, 6)),
    m[, c(2, 5, 6), ])
  expect_equal(
    dust_system_state(obj, index_state = c(1, 2, 4)),
    m[c(1, 2, 4), , ])
  expect_equal(
    dust_system_state(obj, index_group = c(2, 3)),
    m[, , c(2, 3)])
  expect_equal(
    dust_system_state(obj, index_particle = c(4, 7), index_state = c(1, 2, 4)),
    m[c(1, 2, 4), c(4, 7), ])
  expect_equal(
    dust_system_state(obj, index_particle = c(4, 7), index_group = c(2, 3)),
    m[, c(4, 7), c(2, 3)])
  expect_equal(
    dust_system_state(obj, index_state = 2, index_particle = c(4, 7),
                      index_group = c(2, 3)),
    m[2, c(4, 7), c(2, 3), drop = FALSE])

  expect_error(
    dust_system_state(obj, index_particle = integer(0)),
    "'index_particle' must have nonzero length")
  expect_error(
    dust_system_state(obj, index_state = integer(0)),
    "'index_state' must have nonzero length")
  expect_error(
    dust_system_state(obj, index_group = integer(0)),
    "'index_group' must have nonzero length")
})


test_that("can unpack state", {
  set.seed(1)
  pars <- list(len = 5, sd = 1, random_initial = TRUE)
  obj <- dust_system_create(walk(), pars, n_particles = 1)
  dust_system_set_state_initial(obj)
  m <- dust_system_state(obj)
  expect_equal(obj$packer_state$unpack(m),
               list(x = m))
})
