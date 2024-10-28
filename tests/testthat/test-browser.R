test_that("can enable browser", {
  withr::local_options(list(dust.browser_enable = NULL))

  oo <- dust_browser_enabled()
  expect_equal(oo, list("dust.browser_enabled" = NULL))
  expect_true(getOption("dust.browser_enabled"))
  dust_browser_enabled(FALSE)
  expect_false(getOption("dust.browser_enabled"))
  dust_browser_enabled(NULL)
  expect_null(getOption("dust.browser_enabled"))
})


test_that("can change verbosity level", {
  withr::local_options(list(dust.browser_verbosity = NULL))
  oo <- dust_browser_verbosity("verbose")
  expect_equal(oo, list("dust.browser_verbosity" = NULL))
  expect_equal(getOption("dust.browser_verbosity"), "verbose")
  dust_browser_verbosity("quiet")
  expect_equal(getOption("dust.browser_verbosity"), "quiet")
  dust_browser_verbosity("normal")
  expect_equal(getOption("dust.browser_verbosity"), "normal")
  dust_browser_verbosity(NULL)
  expect_null(getOption("dust.browser_verbosity"))
  expect_error(
    dust_browser_verbosity("other"),
    "'level' must be one of 'quiet', 'normal', 'verbose'")
})


test_that("can browser print welcome message", {
  env <- list2env(list(a = 1, b = 2, c = 3))
  withr::with_options(
    list(dust.browser_verbosity = "silent"),
    expect_silent(browser_welcome_message(env)))

  res <- withr::with_options(
    list(dust.browser_verbosity = "normal"),
    evaluate_promise(browser_welcome_message(env, "update", 2)))
  expect_length(res$messages, 1)
  expect_match(
    res$messages,
    "dust debug \\('update'; time = 2\\): see.+dust_browser.+for help")

  res <- withr::with_options(
    list(dust.browser_verbosity = "verbose"),
    evaluate_promise(browser_welcome_message(env, "update", 2)))
  expect_length(res$messages, 5)
  expect_match(
    res$messages[[1]],
    "dust debug \\('update'; time = 2\\)\n")
  expect_match(
    res$messages[[2]],
    "3 variables available")
})


test_that("can't call dust_browser_continue normally", {
  expect_error(
    dust_browser_continue(),
    "Called 'dust_browser_continue()' from outside of a dust browser context",
    fixed = TRUE)
})


test_that("set browser sentinal if requested", {
  skip_if_not_installed("mockery")
  e <- new.env()
  mock_find_env <- mockery::mock(e)
  mockery::stub(dust_browser_continue, "browser_find_parent_env", mock_find_env)
  dust_browser_continue()
  mockery::expect_called(mock_find_env, 1)
  expect_true(e$.dust_browser_continue)
})


test_that("can browse environment", {
  skip_if_not_installed("mockery")
  mock_browse_env <- mockery::mock()
  mockery::stub(browser_env, "browse_env", mock_browse_env)
  env <- new.env()
  env$a <- 1
  withr::with_options(list(dust.browser_enabled = FALSE),
                      browser_env(env, "phase", 1))
  mockery::expect_called(mock_browse_env, 0)

  withr::with_options(
    list(dust.browser_enabled = NULL),
    expect_message(browser_env(env, "phase", 1),
                   "dust debug ('phase'; time = 1):",
                   fixed = TRUE))

  mockery::expect_called(mock_browse_env, 1)
  expect_equal(mockery::mock_args(mock_browse_env)[[1]], list(env))
})
