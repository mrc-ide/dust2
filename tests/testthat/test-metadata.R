test_that("can read sir metadata", {
  meta <- parse_metadata(dust2_file("examples/sir.cpp"))
  expect_equal(meta$class, "sir")
  expect_equal(meta$name, "sir")
  expect_true(meta$has_compare)
  expect_equal(meta$parameters,
               data.frame(name = c("I0", "N", "beta", "gamma", "exp_noise")))
})
