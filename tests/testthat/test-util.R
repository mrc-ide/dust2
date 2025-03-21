test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
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

test_that("fmod works", {
  expect_equal(abs(fmod(5L, 2L)), 1)
  expect_equal(fmod(10, 0.25), 0)
  expect_equal(fmod(10, 0.1), 0)
  expect_equal(fmod(10, 0.01), 0)
})


test_that("describe ranks", {
  expect_equal(rank_description(0), "scalar")
  expect_equal(rank_description(1), "vector")
  expect_equal(rank_description(2), "matrix")
  expect_equal(rank_description(3), "3-dimensional array")
  expect_equal(rank_description(300), "300-dimensional array")
})


test_that("nicely deparse a df", {
  expect_equal(deparse_df(NULL, 2), "NULL")
  expect_equal(deparse_df(data.frame(a = 1), 2),
               'data.frame(\n  a = 1)')
  expect_equal(deparse_df(data.frame(a = 1, b = "2"), 2),
               'data.frame(\n  a = 1,\n  b = "2")')
  expect_equal(deparse_df(data.frame(a = c(1, 2), b = c("3", "4")), 4),
               'data.frame(\n    a = c(1, 2),\n    b = c("3", "4"))')
})
