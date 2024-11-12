test_that("can compile simple system", {
  skip_for_compilation()
  code <- gsub("sir", "mysir", readLines(dust2_file("examples/sir.cpp")))
  filename <- tempfile(fileext = ".cpp")
  writeLines(code, filename)
  gen <- dust_compile(filename, quiet = TRUE, debug = TRUE)
  expect_true(is.function(gen))
  res <- gen()
  expect_s3_class(res, "dust_system_generator")
  expect_equal(attr(res, "name"), "mysir")
  expect_true(attr(res, "properties")$has_compare)

  obj1 <- dust_system_create(gen(), list(), n_particles = 10, seed = 1)
  dust_system_set_state_initial(obj1)
  res1 <- dust_system_simulate(obj1, 0:10)

  obj2 <- dust_system_create(sir(), list(), n_particles = 10, seed = 1)
  dust_system_set_state_initial(obj2)
  expect_equal(res1, dust_system_simulate(obj2, 0:10))

  res <- evaluate_promise(dust_compile(filename, quiet = FALSE, debug = TRUE))
  expect_identical(res$result, gen)
  expect_match(res$messages, "Using cached generator")
})


test_that("can compile into a stable directory", {
  path <- withr::local_tempdir()
  withr::local_envvar(c(DUST_WORKDIR_ROOT = path))

  skip_for_compilation()
  code <- gsub("sir", "mysir2", readLines(dust2_file("examples/sir.cpp")))
  filename <- withr::local_tempfile(fileext = ".cpp")
  writeLines(code, filename)
  gen <- dust_compile(filename, quiet = TRUE, debug = TRUE)
  expect_length(dir(path), 1)

  log <- withr::local_tempfile()
  env <- c(callr::rcmd_safe_env(), DUST_WORKDIR_ROOT = path)
  res <- callr::r(function(filename) {
    dust2::dust_compile(filename, quiet = FALSE, debug = TRUE)
  }, list(filename = filename), env = env, stdout = log, stderr = "2>&1")

  log_txt <- readLines(log)
  expect_match(cli::ansi_strip(log_txt),
               "'src/dust.cpp' is up to date",
               all = FALSE)
  expect_match(cli::ansi_strip(log_txt),
               "Loading mysir2",
               all = FALSE)
  expect_false(any(grepl("compiling", log_txt, ignore.case = TRUE)))
})


test_that("generators can be serialised and used from other processes", {
  skip_if_not_installed("callr")

  code <- gsub("sir", "mysir", readLines(dust2_file("examples/sir.cpp")))
  filename <- tempfile(fileext = ".cpp")
  writeLines(code, filename)
  gen <- dust_compile(filename, quiet = TRUE, debug = TRUE)

  tmp <- withr::local_tempfile()
  saveRDS(gen, tmp)

  log <- withr::local_tempfile()
  expect_equal(
    callr::r(function(path) {
      generator <- readRDS(path)
      dust2::dust_system_state(dust2::dust_system_create(generator, list(), 1))
    }, list(tmp), stdout = log, stderr = "2>&1"),
    numeric(5))
  expect_match(cli::ansi_strip(readLines(log)), "Loading mysir", all = FALSE)
})


test_that("can compile a model with browser-based debugging", {
  expect_no_error(
    dust_compile("examples/debug.cpp", quiet = TRUE, debug = TRUE))
})


test_that("can run system with more variables than output", {
  skip_for_compilation()
  gen <- dust_compile("examples/output.cpp", quiet = TRUE, debug = TRUE)

  sys <- dust_system_create(gen, list())
  dust_system_set_state_initial(sys)
  t <- seq(0, 100, by = 4)
  y <- dust_system_simulate(sys, t)

  sys_cmp <- dust_system_create(sirode, list())
  dust_system_set_state_initial(sys_cmp)
  y_cmp <- dust_system_simulate(sys_cmp, t)

  expect_equal(y[1:4, ], y_cmp[1:4, ], tolerance = 1e-6)
  expect_equal(y[5, ], rep(0, length(t)))
})
