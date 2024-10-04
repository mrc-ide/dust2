test_that("can enable debug", {
  withr::local_options(list(dust.debug_enable = NULL))

  oo <- dust_debug_enabled()
  expect_equal(oo, list("dust.debug_enabled" = NULL))
  expect_true(getOption("dust.debug_enabled"))
  dust_debug_enabled(FALSE)
  expect_false(getOption("dust.debug_enabled"))
  dust_debug_enabled(NULL)
  expect_null(getOption("dust.debug_enabled"))
})


test_that("can change verbosity level", {
  withr::local_options(list(dust.debug_verbosity = NULL))
  oo <- dust_debug_verbosity("verbose")
  expect_equal(oo, list("dust.debug_verbosity" = NULL))
  expect_equal(getOption("dust.debug_verbosity"), "verbose")
  dust_debug_verbosity("quiet")
  expect_equal(getOption("dust.debug_verbosity"), "quiet")
  dust_debug_verbosity("normal")
  expect_equal(getOption("dust.debug_verbosity"), "normal")
  dust_debug_verbosity(NULL)
  expect_null(getOption("dust.debug_verbosity"))
  expect_error(
    dust_debug_verbosity("other"),
    "'level' must be one of 'quiet', 'normal', 'verbose'")
})


test_that("can debug print welcome message", {
  env <- list2env(list(a = 1, b = 2, c = 3))
  withr::with_options(
    list(dust.debug_verbosity = "silent"),
    expect_silent(debug_welcome_message(env)))

  res <- withr::with_options(
    list(dust.debug_verbosity = "normal"),
    evaluate_promise(debug_welcome_message(env)))
  expect_length(res$messages, 1)
  expect_match(
    res$messages,
    "dust debug: see.+dust_debug.+for help")

  res <- withr::with_options(
    list(dust.debug_verbosity = "verbose"),
    evaluate_promise(debug_welcome_message(env)))
  expect_length(res$messages, 4)
  expect_match(
    res$messages[[1]],
    "dust debug\n")
  expect_match(
    res$messages[[2]],
    "3 variables available")
})
