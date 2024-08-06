test_that("Can validate thread count", {
  expect_equal(check_n_threads(4, 4, 1), 4)
  expect_equal(check_n_threads(1, 4, 1), 1)
  expect_equal(check_n_threads(1, 40, 100), 1)
  expect_warning(
    n <- check_n_threads(10, 3, 2),
    "Reducing 'n_threads' from requested 10 to 6")
  expect_equal(n, 6)
})
