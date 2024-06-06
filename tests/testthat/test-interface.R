test_that("can print information about generators", {
  res <- evaluate_promise(withVisible(print(sir())))
  expect_mapequal(res$result, list(value = sir(), visible = FALSE))
  expect_match(res$messages, "<dust_system_generator: sir>",
               fixed = TRUE, all = FALSE)
  expect_match(
    res$messages,
    "Use 'dust2::dust_system_create()' to create a system with this generator",
    fixed = TRUE, all = FALSE)
})


test_that("error if given invalid inputs to dust_system_create", {
  expect_error(
    dust_system_create(NULL),
    "Expected 'generator' to be a 'dust_system_generator' object")
  expect_error(
    dust_system_create("sir"),
    "Expected 'generator' to be a 'dust_system_generator' object")

  foo <- function() {
    dust_system("foo")
  }
  err <- expect_error(
    dust_system_create(foo),
    "Expected 'generator' to be a 'dust_system_generator' object")
  expect_equal(err$body,
               c(i = "Did you mean 'foo()' (i.e., with parentheses)"))
})


test_that("can print information about dust systems", {
  msgs <- function(...) {
    evaluate_promise(print(dust_system_create(sir(), ...)))$messages
  }

  expect_match(msgs(list(), n_particles = 1),
               "5 state x 1 particle\\b",
               all = FALSE)
  expect_match(msgs(list(), n_particles = 10),
               "5 state x 10 particles\\b",
               all = FALSE)
  expect_match(msgs(list(list()), n_particles = 10, n_groups = 1),
               "5 state x 10 particles x 1 group\\b",
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
  expect_error(dust_unfilter_run(NULL),
               "Expected 'unfilter' to be a 'dust_unfilter' object")
  expect_error(dust_filter_run(NULL),
               "Expected 'filter' to be a 'dust_filter' object")
})
