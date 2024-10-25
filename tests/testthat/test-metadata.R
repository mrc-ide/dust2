test_that("can read sir metadata", {
  meta <- parse_metadata(dust2_file("examples/sir.cpp"))
  expect_equal(meta$class, "sir")
  expect_equal(meta$name, "sir")
  expect_true(meta$has_compare)
  expect_equal(meta$default_dt, 1)
  expect_equal(meta$parameters,
               data.frame(name = c("I0", "N", "beta", "gamma", "exp_noise"),
                          type = "real_type",
                          constant = c(FALSE, TRUE, FALSE, FALSE, TRUE),
                          rank = 0L))
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
    "Expected first argument to be an unquoted string",
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


test_that("validate default dt", {
  d <- data.frame(decoration = "dust2::default_dt",
                  params = I(list(alist(0.5))))

  expect_equal(parse_metadata_default_dt(d[-1, ], "discrete"), 1)
  expect_equal(parse_metadata_default_dt(d, "discrete"), 0.5)

  expect_error(
    parse_metadata_default_dt(d, "continuous"),
    "Can't use '[[dust::default_dt()]]' with continuous-time systems",
    fixed = TRUE)

  d$params[[1]] <- list(quote(a))
  expect_error(
    parse_metadata_default_dt(d, "discrete"),
    "Expected a numerical argument to '[[dust2::default_dt()]]'",
    fixed = TRUE)

  d$params[[1]] <- list(dt = 1)
  expect_error(
    parse_metadata_default_dt(d, "discrete"),
    "Expected a single unnamed argument to '[[dust2::default_dt()]]'",
    fixed = TRUE)

  d$params[[1]] <- list(0.37)
  expect_error(
    parse_metadata_default_dt(d, "discrete"),
    "Expected '[[dust2::default_dt()]]' to be the inverse of an integer",
    fixed = TRUE)
})


test_that("validate parameter arguments: name", {
  expect_error(
    parse_metadata_parameter_single(list()),
    "Expected first argument to be an unquoted string")
  expect_error(
    parse_metadata_parameter_single(list(1)),
    "Expected first argument to be an unquoted string")
  expect_error(
    parse_metadata_parameter_single(list(name = 1)),
    "Expected 'name' to be an unquoted string")
  expect_error(
    parse_metadata_parameter_single(list(name = "a")),
    "Expected 'name' to be an unquoted string")
})


test_that("validate parameter arguments: type", {
  expect_error(
    parse_metadata_parameter_single(list(quote(a), type = "foo")),
    "'type' must be one of 'real_type', 'int', 'bool'")
})


test_that("boolean parameters do not have ranges", {
  expect_error(
    parse_metadata_parameter_single(list(quote(a), type = "bool", min = 1)),
    "Argument 'min' does not make sense with 'type = \"bool\"'")
  expect_error(
    parse_metadata_parameter_single(list(quote(a), type = "bool", max = 1)),
    "Argument 'max' does not make sense with 'type = \"bool\"'")
  expect_no_error(
    parse_metadata_parameter_single(list(quote(a), type = "bool", max = NA)))
})


test_that("min, max for numbers must be numbers", {
  expect_error(
    parse_metadata_parameter_single(list(quote(a), min = "1")),
    "Expected 'min' to be numeric")
  expect_error(
    parse_metadata_parameter_single(list(quote(a), max = "1")),
    "Expected 'max' to be numeric")
})


test_that("min and max must leave some range", {
  expect_error(
    parse_metadata_parameter_single(list(quote(a), min = 2, max = 1)),
    "'min' is greater than 'max'")
  expect_no_error(
    parse_metadata_parameter_single(list(quote(a), min = 1, max = 1)))
})


test_that("rank must be a size", {
  expect_no_error(
    parse_metadata_parameter_single(list(quote(a), rank = 0)))
  expect_error(
    parse_metadata_parameter_single(list(quote(a), rank = 1.5)),
    "Expected 'rank' to be integer")
  expect_error(
    parse_metadata_parameter_single(list(quote(a), rank = -1)),
    "'rank' must be at least 0")
})


test_that("required and constant must be booleans", {
  expect_error(
    parse_metadata_parameter_single(list(quote(a), required = "true")),
    "Expected 'required' to be logical")
  expect_error(
    parse_metadata_parameter_single(list(quote(a), constant = "true")),
    "Expected 'constant' to be logical")
})


test_that("if some parameters have required they must all have it", {
  d <- data.frame(line = 1:3,
                  decoration = "dust2::parameter",
                  params = I(list(list(quote(a), required = TRUE),
                                  list(quote(b)),
                                  list(quote(c), required = TRUE))))
  expect_error(
    parse_metadata_parameters(d),
    "Only some '[[dust2::parameter()]]' entries have arguments",
    fixed = TRUE)
})


test_that("prevent duplicate parameter namesx", {
  d <- data.frame(line = 1:3,
                  decoration = "dust2::parameter",
                  params = I(list(list(quote(a)),
                                  list(quote(b)),
                                  list(quote(a)))))
  expect_error(
    parse_metadata_parameters(d),
    "Duplicate parameter names across '[[dust2::paramneter()]]' entries:",
    fixed = TRUE)
})
