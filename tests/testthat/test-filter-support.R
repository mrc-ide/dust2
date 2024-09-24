test_that("can validate that 'dt' is reasonable", {
  expect_no_error(check_dt(1))
  expect_no_error(check_dt(1 / 5))
  expect_error(check_dt(-1, "dt"),
               "Expected 'dt' to be greater than 0")
  expect_error(check_dt(2, "dt"),
               "Expected 'dt' to be at most 1")
  expect_error(check_dt(1 / 3.5, "dt"),
               "Expected 'dt' to be the inverse of an integer")
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
