test_that("can print information about generators", {
  res <- evaluate_promise(withVisible(print(sir())))
  expect_mapequal(res$result, list(value = sir(), visible = FALSE))
  expect_match(res$messages, "<dust_model_generator: sir>",
               fixed = TRUE, all = FALSE)
  expect_match(
    res$messages,
    "Use 'dust2::dust_model_create()' to create a model with this generator",
    fixed = TRUE, all = FALSE)
})


test_that("error if given invalid inputs to dust_model_create", {
  expect_error(
    dust_model_create(NULL),
    "Expected 'generator' to be a 'dust_model_generator' object")
  expect_error(
    dust_model_create("sir"),
    "Expected 'generator' to be a 'dust_model_generator' object")

  foo <- function() {
    dust_model("foo")
  }
  err <- expect_error(
    dust_model_create(foo),
    "Expected 'generator' to be a 'dust_model_generator' object")
  expect_equal(err$body,
               c(i = "Did you mean 'foo()' (i.e., with parentheses)"))
})


test_that("can print information about dust models", {
  msgs <- function(...) {
    evaluate_promise(print(dust_model_create(sir(), ...)))$messages
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
               "This model is deterministic",
               all = FALSE)
})


test_that("error if non-dust model given to dust function", {
  expect_error(dust_model_state(NULL),
               "Expected 'model' to be a 'dust_model' object")
})
