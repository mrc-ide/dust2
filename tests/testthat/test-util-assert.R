test_that("assert_scalar", {
  expect_error(assert_scalar(NULL), "must be a scalar")
  expect_error(assert_scalar(numeric(0)), "must be a scalar")
  expect_error(assert_scalar(1:2), "must be a scalar")
})


test_that("assert_character", {
  expect_silent(assert_character("a"))
  expect_error(assert_character(1), "to be character")
  expect_error(assert_character(TRUE), "to be character")
})


test_that("assert_nonmissing", {
  expect_silent(assert_nonmissing(TRUE))
  expect_error(assert_nonmissing(NA), "Expected 'NA' to be non-NA")
  x <- c(1, NA)
  expect_error(assert_nonmissing(x), "Expected 'x' to be non-NA")
})
