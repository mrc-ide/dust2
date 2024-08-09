test_that("can validate that 'dt' is reasonable", {
  expect_no_error(check_dt(1))
  expect_no_error(check_dt(1 / 5))
  expect_error(check_dt(-1),
               "Expected 'dt' to be greater than 0")
  expect_error(check_dt(2),
               "Expected 'dt' to be at most 1")
  expect_error(check_dt(1 / 3.5),
               "Expected 'dt' to be the inverse of an integer")
})


test_that("can validate time sequence", {
  expect_no_error(check_time_sequence(0, 1:5))

  err <- expect_error(
    check_time_sequence(3, 1:5),
    "Time sequence is not strictly increasing")
  expect_equal(
    err$body,
    c(x = "'time[1]' (1) must be greater than 'time_start' (3)"))

  err <- expect_error(
    check_time_sequence(0, 10:1),
    "Time sequence is not strictly increasing")
  expect_equal(
    err$body,
    c(x = "'time[2]' (9) must be greater than 'time[1]' (10)",
      x = "'time[3]' (8) must be greater than 'time[2]' (9)",
      x = "'time[4]' (7) must be greater than 'time[3]' (8)",
      x = "'time[5]' (6) must be greater than 'time[4]' (7)",
      x = "...and 5 other errors"))

  expect_error(
    check_time_sequence(0, integer()),
    "Expected at least one value in 'time'")
  expect_error(
    check_time_sequence(0, NULL),
    "Expected 'time' to be integer")
})


test_that("check index", {
  expect_no_error(check_index(NULL))
  expect_no_error(check_index(1:4))
  idx <- c(1, 2, 3.4)
  expect_error(check_index(idx),
               "Expected 'idx' to be integer")
  idx <- c(1, 2, -3)
  expect_error(check_index(idx),
               "All elements of 'idx' must be at least 1")
})
