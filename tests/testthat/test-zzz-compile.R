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
