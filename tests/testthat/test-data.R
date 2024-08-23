test_that("can read in simple data set", {
  d <- dust_filter_data(data.frame(time = 1:4, data = 5:8))
  expect_s3_class(d, c("dust_filter_data", "data.frame"), exact = TRUE)
  expect_equal(attr(d, "time"), "time")
  expect_null(attr(d, "group"))
  expect_equal(attr(d, "n_groups"), 1)
})


test_that("can read in simple data set with time in unusual place", {
  expect_error(dust_filter_data(data.frame(t = 1:4, data = 5:8)),
               "Did not find column 'time' in 'data'")
  d <- dust_filter_data(data.frame(t = 1:4, data = 5:8), time = "t")
  expect_equal(attr(d, "time"), "t")
  expect_null(attr(d, "group"))
  expect_equal(attr(d, "n_groups"), 1)
})


test_that("can read in grouped data", {
  data <- cbind(
    expand.grid(time = 1:4, group = 1:3, KEEP.OUT.ATTRS = FALSE),
    data = runif(12))
  d <- dust_filter_data(data)
  expect_s3_class(d, c("dust_fulter_data", "data.frame"))
  expect_equal(attr(d, "time"), "time")
  expect_equal(attr(d, "group"), "group")
  expect_equal(attr(d, "n_groups"), 3)
})


test_that("ignore grouped column for ungrouped data", {
  d <- dust_filter_data(data.frame(time = 1:4, group = "a", data = 5:8))
  expect_equal(attr(d, "time"), "time")
  expect_null(attr(d, "group"))
  expect_equal(attr(d, "n_groups"), 1)
})


test_that("require that data is not empty", {
  expect_error(dust_filter_data(data.frame(time = integer())),
               "Expected 'data' to have at least one row")
  expect_error(
    dust_filter_data(data.frame(time = 1)),
    "Expected 'data' to have at least one column in addition to 'time'")
  expect_error(
    dust_filter_data(data.frame(time = c(1, 1), group = c(1, 2))),
    "Expected 'data' to have at least one column in addition to 'time' and")
})


test_that("time column must be integer-like", {
  expect_error(
    dust_filter_data(data.frame(time = seq(0, 10, length.out = 14))),
    "Expected 'data[[\"time\"]]' to be integer",
    fixed = TRUE)
})


test_that("require group if times duplicate", {
  d <- data.frame(time = c(1, 1, 2, 2), y = 1:4)
  expect_error(
    dust_filter_data(d),
    "'data' looks grouped, but 'group' not specified")
  expect_error(
    dust_filter_data(d, group = "group"),
    "Did not find column 'group' in 'data'")
  expect_error(
    dust_filter_data(d, group = "other"),
    "Did not find column 'other' in 'data'")
})


test_that("validate that group is usable", {
  d <- data.frame(time = c(1, 1, 2, 2),
                  nm = c("a", "a", "b", "b"),
                  skip = c(1, 3, 1, 3),
                  imbalanced = c(1, 2, 1, 1))
  expect_error(
    dust_filter_data(d, group = "nm"),
    "Expected 'data[[\"nm\"]]' to be integer",
    fixed = TRUE)
  expect_error(
    dust_filter_data(d, group = "skip"),
    "Expected 'data[[\"skip\"]]' to contain integer values in [1, 2]",
    fixed = TRUE)
  expect_error(
    dust_filter_data(d, group = "imbalanced"),
    "Not all groups in 'data' have the same times")
})


test_that("if an object has class 'dust_filter_data' we just pass it back", {
  x <- structure(list(), class = "dust_filter_data")
  expect_identical(dust_filter_data(x), x)
})


test_that("can convert a data.frame with list columns", {
  y <- lapply(1:5, function(i) runif(2))
  d <- data.frame(time = 1:5, x = 1:5, y = I(y),
                  stringsAsFactors = FALSE)
  d2 <- prepare_data(d, 1)
  expect_equal(d2$data[[1]], list(list(x = 1, y = y[[1]])))
  expect_equal(d2$data[[3]], list(list(x = 3, y = y[[3]])))
})
