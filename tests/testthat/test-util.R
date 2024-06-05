test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})


test_that("can error if packages are not installed", {
  skip_if_not_installed("mockery")
  mock_require_ns <- mockery::mock(TRUE, FALSE, FALSE)
  mockery::stub(stop_unless_installed, "requireNamespace", mock_require_ns)
  expect_error(
    stop_unless_installed(c("a", "b", "c")),
    "Packages missing for this functionality: 'b' and 'c'")
  mockery::expect_called(mock_require_ns, 3)
})
