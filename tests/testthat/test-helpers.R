test_that("validate dimensions of a vector", {
  expect_true(test_check_dimensions(1:5, 5L, "x"))
  expect_error(
    test_check_dimensions(1:5, 6L, "x"),
    "Expected 'x' to have length 6, but it had length 5")
  expect_error(
    test_check_dimensions(cbind(1:5), 5L, "x"),
    "Expected 'x' to be a vector, but was given a matrix")
  expect_error(
    test_check_dimensions(array(1:5, c(1, 1, 5)), 5L, "x"),
    "Expected 'x' to be a vector, but was given a 3-dimensional array")
})


test_that("validate dimensions of a matrix", {
  expect_true(test_check_dimensions(matrix(0, 2, 3), c(2L, 3L), "x"))
  expect_error(
    test_check_dimensions(matrix(0, 2, 3), c(5L, 3L), "x"),
    "Expected 'x' to have 5 rows, but was given 2")
  expect_error(
    test_check_dimensions(matrix(0, 2, 3), c(2L, 6L), "x"),
    "Expected 'x' to have 6 columns, but was given 3")
  expect_error(
    test_check_dimensions(1:6, c(2L, 3L), "x"),
    "Expected 'x' to be a matrix, but was given a vector")
  expect_error(
    test_check_dimensions(array(1:6, c(1, 2, 3)), c(2L, 3L), "x"),
    "Expected 'x' to be a matrix, but was given a 3-dimensional array")
})


test_that("validate dimensions of an array", {
  expect_true(test_check_dimensions(array(0, 2:4), c(2L, 3L, 4L), "x"))
  expect_error(
    test_check_dimensions(array(0, 2:4), c(5L, 3L, 4L), "x"),
    "Expected dimension 1 of 'x' to have length 5, but it had length 2")
  expect_error(
    test_check_dimensions(array(0, 2:4), c(2L, 6L, 4L), "x"),
    "Expected dimension 2 of 'x' to have length 6, but it had length 3")
  expect_error(
    test_check_dimensions(1:24, 2:4, "x"),
    "Expected 'x' to be a 3-dimensional array, but was given a vector")
  expect_error(
    test_check_dimensions(matrix(0, c(3, 8)), 2:4, "x"),
    "Expected 'x' to be a 3-dimensional array, but was given a matrix")
  expect_error(
    test_check_dimensions(array(0, 1:4), 2:4, "x"),
    "Expected 'x' to be a 3-dimensional array, but was given a 4-dimensional")
})


test_that("can read dimensions from vector", {
  expect_equal(test_read_dimensions(0:3, 1, "x"), 4)
  expect_error(test_read_dimensions(cbind(0:3), 1, "x"),
               "Expected 'x' to be a vector, but was given a matrix")
  expect_error(
    test_read_dimensions(array(0:3, c(1, 2, 2)), 1, "x"),
    "Expected 'x' to be a vector, but was given a 3-dimensional array")
})


test_that("can read dimensions from vector", {
  expect_equal(test_read_dimensions(matrix(0, 2, 3), 2, "x"), c(2, 3))
  expect_error(test_read_dimensions(1:6, 2, "x"),
               "Expected 'x' to be a matrix, but was given a vector")
  expect_error(
    test_read_dimensions(array(0, 1:3), 2, "x"),
    "Expected 'x' to be a matrix, but was given a 3-dimensional array")
})


test_that("can read dimensions from array", {
  expect_equal(test_read_dimensions(array(0, 2:4), 3, "x"), 2:4)
  expect_error(
    test_read_dimensions(1:6, 3, "x"),
    "Expected 'x' to be a 3-dimensional array, but was given a vector")
  expect_error(
    test_read_dimensions(matrix(0, 2, 3), 3, "x"),
    "Expected 'x' to be a 3-dimensional array, but was given a matrix")
  expect_error(
    test_read_dimensions(array(0, 1:4), 3, "x"),
    "Expected 'x' to be a 3-dimensional array, but was given a 4-dimensional")
})
