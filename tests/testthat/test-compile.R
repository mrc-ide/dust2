test_that("can compile simple model", {
  code <- gsub("sir", "mysir", readLines(dust2_file("examples/sir.cpp")))
  filename <- tempfile(fileext = ".cpp")
  writeLines(code, filename)
  gen <- dust_compile(filename, quiet = TRUE, debug = TRUE)
  expect_true(is.function(gen))
  res <- gen()
  expect_s3_class(res, "dust_model_generator")
  expect_equal(res$name, "mysir")
  expect_true(res$properties$has_compare)

  obj1 <- dust_model_create(gen(), list(), n_particles = 10, seed = 1)
  dust_model_set_state_initial(obj1)
  res1 <- dust_model_simulate(obj1, 0:10)

  obj2 <- dust_model_create(sir(), list(), n_particles = 10, seed = 1)
  dust_model_set_state_initial(obj2)
  expect_equal(res1, dust_model_simulate(obj2, 0:10))
})
