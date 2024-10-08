test_that("binary search", {
  t <- 0:9 + 0.5
  eps <- 1e-7

  expect_equal(test_interpolate_search(0.5, t), 0)
  expect_equal(test_interpolate_search(0.5 + eps, t), 0)
  expect_equal(test_interpolate_search(1.5 - eps, t), 0)

  expect_equal(test_interpolate_search(1.5, t), 1)
  expect_equal(test_interpolate_search(1.5 + eps, t), 1)
  expect_equal(test_interpolate_search(2.5 - eps, t), 1)

  expect_equal(test_interpolate_search(5.5, t), 5)
  expect_equal(test_interpolate_search(5.5 + eps, t), 5)
  expect_equal(test_interpolate_search(6.5 - eps, t), 5)

  expect_equal(test_interpolate_search(9.5 - eps, t), 8)
  expect_equal(test_interpolate_search(9.5, t), 9)

  expect_error(
    test_interpolate_search(0, t),
    paste("Tried to interpolate at time = 0, which is 0.5",
          "before the first time (0.5)"),
    fixed = TRUE)
  expect_error(
    test_interpolate_search(9.5 + eps, t),
    paste("Tried to interpolate at time = .+",
          "which is .+ after the last time \\(9.5\\)"))
})


test_that("can work with simple constant interpolation", {
  set.seed(1)
  t <- as.numeric(0:10)
  y <- runif(length(t))
  expect_error(
    test_interpolate_constant1(t, y, 0 - 1e-8),
    "Tried to interpolate.+before the first time")
  expect_equal(test_interpolate_constant1(t, y, 0), y[[1]])
  expect_equal(test_interpolate_constant1(t, y, 1 - 1e-8), y[[1]])
  expect_equal(test_interpolate_constant1(t, y, 1), y[[2]])
  expect_equal(test_interpolate_constant1(t, y, 2), y[[3]])

  z <- vapply(t, function(z) test_interpolate_constant1(t, y, z), numeric(1))
  expect_equal(z, y)

  expect_equal(test_interpolate_constant1(t, y, 10 - 1e-8), y[[10]])
  expect_equal(test_interpolate_constant1(t, y, 10), y[[11]])
  expect_equal(test_interpolate_constant1(t, y, 100), y[[11]])
})


test_that("can work with simple linear interpolation", {
  set.seed(1)
  t <- as.numeric(0:10)
  y <- runif(length(t))
  expect_error(
    test_interpolate_linear1(t, y, 0 - 1e-8),
    "Tried to interpolate.+before the first time")
  expect_error(
    test_interpolate_linear1(t, y, 10 + 1e-8),
    "Tried to interpolate.+after the last time")
  cmp <- approxfun(t, y)

  expect_equal(test_interpolate_linear1(t, y, 0), y[[1]])
  expect_equal(test_interpolate_linear1(t, y, 1 - 1e-8), cmp(1 - 1e-8))
  expect_equal(test_interpolate_linear1(t, y, 0.5), cmp(0.5))
  expect_equal(test_interpolate_linear1(t, y, 1), y[[2]])
  expect_equal(test_interpolate_linear1(t, y, 2), y[[3]])

  z <- vapply(t, function(z) test_interpolate_linear1(t, y, z), numeric(1))
  expect_equal(z, y)
})


test_that("can work with simple spline interpolation", {
  set.seed(1)
  t <- as.numeric(0:10)
  y <- runif(length(t))

  expect_error(
    test_interpolate_spline1(t, y, 0 - 1e-8),
    "Tried to interpolate.+before the first time")
  expect_error(
    test_interpolate_spline1(t, y, 10 + 1e-8),
    "Tried to interpolate.+after the last time")
  cmp <- splinefun(t, y, method = "natural")

  z <- vapply(t, function(z) test_interpolate_spline1(t, y, z), numeric(1))
  expect_equal(z, y)

  expect_equal(test_interpolate_spline1(t, y, 0), y[[1]])
  expect_equal(test_interpolate_spline1(t, y, 1 - 1e-8), cmp(1 - 1e-8))
  expect_equal(test_interpolate_spline1(t, y, 0.5), cmp(0.5))
  expect_equal(test_interpolate_spline1(t, y, 1), y[[2]])
  expect_equal(test_interpolate_spline1(t, y, 2), y[[3]])
})


test_that("Check that time values are sensible", {
  t <- c(0, 1, 1, 2)
  y <- c(0, 1, 2, 3)
  expect_error(
    test_interpolate_spline1(t, y, 1),
    "Time variable 't' must be strictly increasing but was not at index 2")
  expect_error(
    test_interpolate_spline1(t[-2], y, 1),
    "Time variable 't' and interpolation target 'y' must have the same length")
})


test_that("can interpolate a matrix, piecewise constant", {
  t <- c(0, 1, 2, 3, 4)
  nt <- length(t)
  nr <- 7
  nc <- 3
  y <- array(runif(nt * nr * nc), c(nr, nc, nt))
  expect_equal(test_interpolate_constant2(t, y, 0), c(y[, , 1]))
  expect_equal(test_interpolate_constant2(t, y, 1 - 1e-6), c(y[, , 1]))
  expect_equal(test_interpolate_constant2(t, y, 1), c(y[, , 2]))
  expect_equal(test_interpolate_constant2(t, y, 10), c(y[, , 5]))
})


test_that("can interpolate a matrix, piecewise linear", {
  t <- c(0, 1, 2, 3, 4)
  nt <- length(t)
  nr <- 7
  nc <- 3
  y <- array(runif(nt * nr * nc), c(nr, nc, nt))
  expect_equal(test_interpolate_linear2(t, y, 0), c(y[, , 1]))
  expect_equal(test_interpolate_linear2(t, y, 1), c(y[, , 2]))
  expect_equal(test_interpolate_linear2(t, y, 4), c(y[, , 5]))
  expect_equal(test_interpolate_linear2(t, y, 1.3),
               c(y[, , 2] + (y[, , 3] - y[, , 2]) * 0.3))
})


test_that("can interpolate a matrix, spline", {
  set.seed(1)
  t <- c(0, 1, 2, 3, 4)
  nt <- length(t)
  nr <- 7
  nc <- 3
  y <- array(runif(nt * nr * nc), c(nr, nc, nt))

  expect_equal(test_interpolate_spline2(t, y, 0), c(y[, , 1]))
  expect_equal(test_interpolate_spline2(t, y, 1), c(y[, , 2]))
  expect_equal(test_interpolate_spline2(t, y, 4), c(y[, , 5]))

  yy <- apply(array(y, c(nr * nc, nt)), 1,
              function(yi) test_interpolate_spline1(t, yi, 1.3))
  expect_equal(test_interpolate_spline2(t, y, 1.3),
               yy)
})
