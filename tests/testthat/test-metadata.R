test_that("can read sir metadata", {
  meta <- parse_metadata(dust2_file("examples/sir.cpp"))
  expect_equal(meta$class, "sir")
  expect_equal(meta$name, "sir")
  expect_true(meta$has_compare)
  expect_equal(meta$default_dt, 1)
  expect_equal(meta$parameters,
               data.frame(name = c("I0", "N", "beta", "gamma", "exp_noise")))
})


test_that("can validate class metadata", {
  tmp <- withr::local_tempfile()
  writeLines("// [[dust2::class(a, b)]]", tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected a single unnamed argument to '[[dust2::class()]]'",
    fixed = TRUE)

  writeLines("// [[dust2::class('x')]]", tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected an unquoted string argument to '[[dust2::class()]]'",
    fixed = TRUE)

  writeLines("// [[dust2::klass(a, b)]]", tmp)
  expect_error(
    parse_metadata(tmp),
    "Attribute '[[dust2::class()]]' is required, but was not found",
    fixed = TRUE)

  writeLines(c("// [[dust2::class(x)]]", "// [[dust2::class(y)]]"), tmp)
  expect_error(
    parse_metadata(tmp),
    "More than one '[[dust2::class()]]' attribute found",
    fixed = TRUE)

  writeLines(c("// [[dust2::class(a)]]", "// [[dust2::time_type(discrete)]]"),
             tmp)
  expect_equal(parse_metadata(tmp)$class, "a")
})


test_that("can validate name metadata", {
  base <- c("// [[dust2::class(a)]]", "// [[dust2::time_type(discrete)]]")
  tmp <- withr::local_tempfile()
  writeLines(c(base, "// [[dust2::name()]]"),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected a single unnamed argument to '[[dust2::name()]]'",
    fixed = TRUE)

  writeLines(c(base, "// [[dust2::name('b')]]"),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected an unquoted string argument to '[[dust2::name()]]'",
    fixed = TRUE)

  writeLines(base, tmp)
  expect_equal(parse_metadata(tmp)$name, "a")
  writeLines(c(base, "// [[dust2::name(b)]]"), tmp)
  expect_equal(parse_metadata(tmp)$name, "b")
})


test_that("can validate time type metadata", {
  base <- "// [[dust2::class(a)]]"
  tmp <- withr::local_tempfile()

  writeLines(c(base, "// [[dust2::time_type()]]"),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected a single unnamed argument to '[[dust2::time_type()]]'",
    fixed = TRUE)

  writeLines(c(base, '// [[dust2::time_type("discrete")]]'),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected an unquoted string argument to '[[dust2::time_type()]]'",
    fixed = TRUE)

  writeLines(c(base, "// [[dust2::time_type(bouncy)]]"),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected argument to '[[dust2::time_type()]]' to be one of 'discrete'",
    fixed = TRUE)

  writeLines(c(base, "// [[dust2::time_type(discrete)]]"),
             tmp)
  expect_equal(parse_metadata(tmp)$time_type, "discrete")
})


test_that("can validate compare metadata", {
  base <- c("// [[dust2::class(a)]]", "// [[dust2::time_type(discrete)]]")
  tmp <- withr::local_tempfile()
  writeLines(c(base, "// [[dust2::has_compare(TRUE)]]"),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected no arguments to '[[dust2::has_compare()]]'",
    fixed = TRUE)

  writeLines(base, tmp)
  expect_false(parse_metadata(tmp)$has_compare)
  writeLines(c(base, "// [[dust2::has_compare()]]"), tmp)
  expect_true(parse_metadata(tmp)$has_compare)
})


test_that("can validate parameter metdata", {
  base <- c("// [[dust2::class(a)]]", "// [[dust2::time_type(discrete)]]")
  tmp <- withr::local_tempfile()
  writeLines(c(base, "// [[dust2::parameter(TRUE)]]"),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected an unnamed unquoted string argument to",
    fixed = TRUE)
})


test_that("require that file exists", {
  expect_error(
    parse_metadata(tempfile()),
    "File '.+' does not exist")
})


test_that("can specify default dt in discrete time models", {
  tmp <- withr::local_tempfile()
  writeLines(c(
    "// [[dust2::class(a)]]",
    "// [[dust2::time_type(discrete)]]",
    "// [[dust2::default_dt(0.25)]]"),
    tmp)
  expect_equal(parse_metadata(tmp)$default_dt, 0.25)
})


test_that("can validate default dt in discrete time models", {
  tmp <- withr::local_tempfile()
  writeLines(c(
    "// [[dust2::class(a)]]",
    "// [[dust2::time_type(discrete)]]",
    "// [[dust2::default_dt(0.32)]]"),
    tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected '[[dust2::default_dt()]]' to be the inverse of an integer",
    fixed = TRUE)
})
