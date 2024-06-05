test_that("can read sir metadata", {
  meta <- parse_metadata(dust2_file("examples/sir.cpp"))
  expect_equal(meta$class, "sir")
  expect_equal(meta$name, "sir")
  expect_true(meta$has_compare)
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

  writeLines("// [[dust2::class(a)]]", tmp)
  expect_equal(parse_metadata(tmp)$class, "a")
})


test_that("can validate name metadata", {
  tmp <- withr::local_tempfile()
  writeLines(c("// [[dust2::class(a)]]", "// [[dust2::name()]]"),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected a single unnamed argument to '[[dust2::name()]]'",
    fixed = TRUE)

  writeLines(c("// [[dust2::class(a)]]", "// [[dust2::name('b')]]"),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected an unquoted string argument to '[[dust2::name()]]'",
    fixed = TRUE)

  writeLines("// [[dust2::class(a)]]", tmp)
  expect_equal(parse_metadata(tmp)$name, "a")
  writeLines(c("// [[dust2::class(a)]]", "// [[dust2::name(b)]]"), tmp)
  expect_equal(parse_metadata(tmp)$name, "b")
})


test_that("can validate compare metadata", {
  tmp <- withr::local_tempfile()
  writeLines(c("// [[dust2::class(a)]]", "// [[dust2::has_compare(TRUE)]]"),
             tmp)
  expect_error(
    parse_metadata(tmp),
    "Expected no arguments to '[[dust2::has_compare()]]'",
    fixed = TRUE)

  writeLines("// [[dust2::class(a)]]", tmp)
  expect_false(parse_metadata(tmp)$has_compare)
  writeLines(c("// [[dust2::class(a)]]", "// [[dust2::has_compare()]]"), tmp)
  expect_true(parse_metadata(tmp)$has_compare)
})


test_that("can validate parameter metdata", {
  tmp <- withr::local_tempfile()
  writeLines(c("// [[dust2::class(a)]]", "// [[dust2::parameter(TRUE)]]"),
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
