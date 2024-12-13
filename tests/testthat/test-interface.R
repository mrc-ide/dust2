test_that("can print information about generators", {
  res <- evaluate_promise(withVisible(print(sir())))
  expect_mapequal(res$result, list(value = sir(), visible = FALSE))
  expect_match(res$messages, "<dust_system_generator: sir>",
               fixed = TRUE, all = FALSE)
  expect_match(
    res$messages,
    "Use .*dust2::dust_system_create.* to create a system with this generator",
    all = FALSE)
})


test_that("error if given invalid inputs to dust_system_create", {
  expect_error(
    dust_system_create(NULL),
    "Expected 'generator' to be a 'dust_system_generator' object")
  expect_error(
    dust_system_create("sir"),
    "Expected 'generator' to be a 'dust_system_generator' object")
})


test_that("can print information about dust systems", {
  msgs <- function(...) {
    evaluate_promise(print(dust_system_create(sir(), ...)))$messages
  }

  expect_match(msgs(list(), n_particles = 1),
               "single particle with 5 states\\b",
               all = FALSE)
  expect_match(msgs(list(), n_particles = 1,
                    preserve_particle_dimension = TRUE),
               "5 state x 1 particle\\b",
               all = FALSE)
  expect_match(msgs(list(), n_particles = 10),
               "5 state x 10 particles\\b",
               all = FALSE)
  expect_match(msgs(list(list()), n_particles = 10, n_groups = 1,
                    preserve_group_dimension = TRUE),
               "5 state x 10 particles x 1 group\\b",
               all = FALSE)
  expect_match(msgs(list(list()), n_particles = 10, n_groups = 1),
               "5 state x 10 particles\\b",
               all = FALSE)
  expect_match(msgs(list(list(), list()), n_particles = 10, n_groups = 2),
               "5 state x 10 particles x 2 groups\\b",
               all = FALSE)
  expect_match(msgs(list(), n_particles = 1, deterministic = TRUE),
               "This system is deterministic",
               all = FALSE)
})


test_that("error if non-dust system given to dust function", {
  expect_error(dust_system_state(NULL),
               "Expected 'sys' to be a 'dust_system' object")
})


test_that("error if non-dust system given to dust function", {
  thing <- NULL
  expect_error(dust_likelihood_run(thing),
               "Expected 'obj' to be a 'dust_likelihood' object")
})


test_that("can control dropping of single-element group dimensions", {
  pars <- list(len = 3, sd = 1, random_initial = TRUE)
  obj1 <- dust_system_create(walk(), list(pars), n_particles = 10,
                             preserve_group_dimension = TRUE,
                             seed = 42)
  obj2 <- dust_system_create(walk(), pars, n_particles = 10,
                             seed = 42)
  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)
  rng <- monty::monty_rng_create(seed = 42, n_streams = 10)
  r <- monty::monty_random_n_normal(3, 0, 1, rng)
  expect_equal(dim(obj1), c(3, 10, 1))
  expect_equal(dim(obj2), c(3, 10))

  expect_equal(dust_system_state(obj1), array(r, c(3, 10, 1)))
  expect_equal(dust_system_state(obj2), r)
})


test_that("can control dropping of single-element particle dimensions", {
  pars <- list(list(len = 3, sd = 1, random_initial = TRUE),
               list(len = 3, sd = 1, random_initial = TRUE))
  obj1 <- dust_system_create(walk(), pars, n_particles = 1, n_groups = 2,
                             preserve_particle_dimension = TRUE,
                             seed = 42)
  obj2 <- dust_system_create(walk(), pars, n_particles = 1, n_groups = 2,
                             seed = 42)
  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)
  rng <- monty::monty_rng_create(seed = 42, n_streams = 2)
  r <- monty::monty_random_n_normal(3, 0, 1, rng)
  expect_equal(dim(obj1), c(3, 1, 2))
  expect_equal(dim(obj2), c(3, 2))

  expect_equal(dust_system_state(obj1), array(r, c(3, 1, 2)))
  expect_equal(dust_system_state(obj2), r)
})


test_that("can control dropping of all dimensions", {
  pars <- list(len = 3, sd = 1, random_initial = TRUE)
  obj1 <- dust_system_create(walk(), list(pars), n_particles = 1, n_groups = 1,
                             preserve_group_dimension = TRUE,
                             preserve_particle_dimension = TRUE,
                             seed = 42)
  obj2 <- dust_system_create(walk(), pars, n_particles = 1, n_groups = 1,
                             seed = 42)

  dust_system_set_state_initial(obj1)
  dust_system_set_state_initial(obj2)
  rng <- monty::monty_rng_create(seed = 42, n_streams = 1)
  r <- monty::monty_random_n_normal(3, 0, 1, rng)
  expect_equal(dim(obj1), c(3, 1, 1))
  expect_equal(dim(obj2), 3)

  expect_equal(dust_system_state(obj1), array(r, c(3, 1, 1)))
  expect_equal(dust_system_state(obj2), drop(r))
})


test_that("format dimensions as a string", {
  expect_equal(
    format_dimensions(list(preserve_particle_dimension = TRUE,
                           preserve_group_dimension = TRUE,
                           n_particles = 5,
                           n_groups = 3,
                           n_state = 10)),
    "10 state x 5 particles x 3 groups")
  expect_equal(
    format_dimensions(list(preserve_particle_dimension = TRUE,
                           preserve_group_dimension = TRUE,
                           n_particles = 5,
                           n_groups = 3)),
    "5 particles x 3 groups")

  expect_equal(
    format_dimensions(list(preserve_particle_dimension = FALSE,
                           preserve_group_dimension = TRUE,
                           n_particles = 1,
                           n_groups = 3,
                           n_state = 10)),
    "10 state x 3 groups")
  expect_equal(
    format_dimensions(list(preserve_particle_dimension = FALSE,
                           preserve_group_dimension = TRUE,
                           n_particles = 1,
                           n_groups = 3)),
    "3 groups")

  expect_equal(
    format_dimensions(list(preserve_particle_dimension = TRUE,
                           preserve_group_dimension = FALSE,
                           n_particles = 5,
                           n_groups = 1,
                           n_state = 10)),
    "10 state x 5 particles")
  expect_equal(
    format_dimensions(list(preserve_group_dimension = FALSE,
                           n_particles = 5,
                           n_groups = 1,
                           n_state = 10)),
    "10 state x 5 particles")
  expect_equal(
    format_dimensions(list(preserve_particle_dimension = TRUE,
                           preserve_group_dimension = FALSE,
                           n_particles = 5,
                           n_groups = 1)),
    "5 particles")

  expect_equal(
    format_dimensions(list(preserve_particle_dimension = FALSE,
                           preserve_group_dimension = FALSE,
                           n_particles = 1,
                           n_groups = 1,
                           n_state = 10)),
    "single particle with 10 states")
  expect_equal(
    format_dimensions(list(preserve_particle_dimension = FALSE,
                           preserve_group_dimension = FALSE,
                           n_particles = 1,
                           n_groups = 1)),
    "single particle")
})


test_that("can describe time", {
  expect_equal(
    describe_time("continuous", NULL),
    "This system runs in continuous time")
  expect_equal(
    describe_time("discrete", 1),
    "This system runs in discrete time with a default dt of 1")
  expect_equal(
    describe_time("mixed", NULL),
    paste("This system runs in both continuous time and discrete time",
          "with discrete time disabled by default"))
  expect_equal(
    describe_time("mixed", 1),
    paste("This system runs in both continuous time and discrete time",
          "with a default dt of 1"))

  expect_equal(
    describe_time("continuous", NULL, NULL),
    "This system runs in continuous time")
  expect_equal(
    describe_time("discrete", 1, 0.5),
    "This system runs in discrete time with dt = 0.5")
  expect_equal(
    describe_time("mixed", NULL, NULL),
    paste("This system runs in both continuous time and discrete time",
          "with discrete time disabled"))
  expect_equal(
    describe_time("mixed", NULL, 0.5),
    paste("This system runs in both continuous time and discrete time",
          "with dt = 0.5"))
})


test_that("system objects are immutable", {
  obj <- dust_system_create(walk, list(sd = 1))
  expect_error(
    obj$n_threads <- 1,
    "Cannot write to 'dust_system' objects, they are read-only")
  expect_error(
    obj[["n_threads"]] <- 1,
    "Cannot write to 'dust_system' objects, they are read-only")
  expect_error(
    obj[[1]] <- 1,
    "Cannot write to 'dust_system' objects, they are read-only")
  expect_error(
    obj[1] <- 1,
    "Cannot write to 'dust_system' objects, they are read-only")
})


test_that("system generator objects are immutable", {
  expect_error(
    sir$n_threads <- 1,
    "Cannot write to 'dust_system_generator' objects, they are read-only")
  expect_error(
    sir[["n_threads"]] <- 1,
    "Cannot write to 'dust_system_generator' objects, they are read-only")
  expect_error(
    sir[[1]] <- 1,
    "Cannot write to 'dust_system_generator' objects, they are read-only")
  expect_error(
    sir[1] <- 1,
    "Cannot write to 'dust_system_generator' objects, they are read-only")
})
