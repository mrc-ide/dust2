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


test_that("writelines_if_changed doesn't replace file", {
  workdir <- withr::local_tempdir()
  path <- "myfile"
  path_full <- file.path(workdir, path)
  text1 <- c("a", "b", "c")
  text2 <- c("a", "b", "c", "d")
  str1 <- paste(text1, collapse = "\n")
  str2 <- structure(str1, class = "glue")
  expect_silent(writelines_if_changed(text1, workdir, path, quiet = TRUE))
  expect_true(same_content(path_full, text1))
  expect_true(same_content(path_full, str1))
  expect_true(same_content(path_full, str2))
  expect_silent(writelines_if_changed(str1, workdir, path, quiet = TRUE))
  expect_silent(writelines_if_changed(str2, workdir, path, quiet = TRUE))
  expect_true(file.exists(path_full))
  expect_equal(readLines(path_full), text1)
  t <- file.mtime(path_full)
  expect_silent(writelines_if_changed(text1, workdir, path, quiet = TRUE))
  expect_identical(file.mtime(path_full), t)
  expect_silent(writelines_if_changed(text2, workdir, path, quiet = TRUE))
  expect_equal(readLines(path_full), text2)
  ## I don't trust times and sub-second accuracy not guaranted; see
  ## ?file.mtime
  skip_on_cran()
  skip_on_os("windows")
  expect_gt(file.mtime(path_full), t)
})


test_that("tell user about changes to files", {
  workdir <- withr::local_tempdir()

  text1 <- c("a", "b", "c")
  text2 <- c("a", "b", "c", "d")
  str1 <- paste(text1, collapse = "\n")
  str2 <- paste(text2, collapse = "\n")
  expect_message(
    writelines_if_changed(text1, workdir, "myfile", FALSE),
    "Wrote 'myfile'")
  expect_message(
    writelines_if_changed(text1, workdir, "myfile", FALSE),
    "'myfile' is up to date")
  expect_message(
    writelines_if_changed(text2, workdir, "myfile", FALSE),
    "Wrote 'myfile'")
})


test_that("can set names onto things", {
  expect_null(set_names(NULL, "x"))
  expect_equal(set_names("a", "x"), c(x = "a"))
  expect_equal(set_names(c("a", "b"), "x"), c(x = "a", x = "b"))
  expect_equal(set_names(c("a", "b"), c("x", "i")), c(x = "a", i = "b"))
})


test_that("can use protect to avoid errors", {
  f <- function(x) {
    if (x < 0) {
      stop("invalid value for 'x'")
    }
    sqrt(x)
  }
  g <- protect(f, -Inf)
  expect_equal(g(4), 2)
  expect_equal(g(-1), -Inf)
})
