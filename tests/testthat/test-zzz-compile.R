test_that("can compile simple system", {
  skip_for_compilation()
  code <- gsub("sir", "mysir", readLines(dust2_file("examples/sir.cpp")))
  filename <- tempfile(fileext = ".cpp")
  writeLines(code, filename)
  gen <- dust_compile(filename, quiet = TRUE, debug = TRUE)
  expect_true(is.function(gen))
  res <- gen()
  expect_s3_class(res, "dust_system_generator")
  expect_equal(res$name, "mysir")
  expect_true(res$properties$has_compare)

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
  expect_match(log_txt, "'src/dust.cpp' is up to date", all = FALSE)
  expect_match(log_txt, "Loading(.+)mysir2", all = FALSE)
  expect_false(any(grepl("compiling", log_txt, ignore.case = TRUE)))
})
