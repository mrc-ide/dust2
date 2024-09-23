test_that("can sum over vectors", {
  x <- runif(10)
  expect_equal(test_sum(x, NULL), sum(x))
  expect_equal(test_sum(x, c(4L, 6L)), sum(x[4:6]))
  expect_equal(test_sum(x, c(5L, 5L)), x[5])
  expect_equal(test_sum(x, c(6L, 4L)), 0)
})


test_that("can sum over matrices", {
  x <- matrix(runif(15), 3, 5)
  expect_equal(test_sum(x, NULL), sum(x))
  expect_equal(test_sum(x, c(1L, 2L, 3L, 4L)), sum(x[1:2, 3:4]))
  expect_equal(test_sum(x, c(2L, 2L, 2L, 4L)), sum(x[2, 2:4]))
  expect_equal(test_sum(x, c(2L, 2L, 4L, 2L)), 0)
})


test_that("can take min over 3d arrays", {
  x <- array(runif(3 * 5 * 7), c(3, 5, 7))
  expect_equal(test_min(x, NULL), min(x))
  expect_equal(test_min(x, c(1L, 2L, 3L, 4L, 5L, 6L)),
               min(x[1:2, 3:4, 5:6]))
  expect_equal(test_min(x, c(2L, 2L, 2L, 4L, 1L, 1L)), min(x[2, 2:4, 1]))
  expect_equal(test_min(x, c(2L, 2L, 4L, 2L, 1L, 1L)), 0)
})


test_that("can take min over vectors", {
  x <- runif(10)
  expect_equal(test_min(x, NULL), min(x))
  expect_equal(test_min(x, c(4L, 6L)), min(x[4:6]))
  expect_equal(test_min(x, c(5L, 5L)), x[5])
  expect_equal(test_min(x, c(6L, 4L)), Inf)
})


test_that("can take min over matrices", {
  x <- matrix(runif(15), 3, 5)
  expect_equal(test_min(x, NULL), min(x))
  expect_equal(test_min(x, c(1L, 2L, 3L, 4L)), min(x[1:2, 3:4]))
  expect_equal(test_min(x, c(2L, 2L, 2L, 4L)), min(x[2, 2:4]))
  expect_equal(test_min(x, c(2L, 2L, 4L, 2L)), Inf)
})


test_that("can take min over 3d arrays", {
  x <- array(runif(3 * 5 * 7), c(3, 5, 7))
  expect_equal(test_min(x, NULL), min(x))
  expect_equal(test_min(x, c(1L, 2L, 3L, 4L, 5L, 6L)),
               min(x[1:2, 3:4, 5:6]))
  expect_equal(test_min(x, c(2L, 2L, 2L, 4L, 1L, 1L)), min(x[2, 2:4, 1]))
  expect_equal(test_min(x, c(2L, 2L, 4L, 2L, 1L, 1L)), Inf)
})
