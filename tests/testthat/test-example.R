test_that("Can create example", {
  expect_identical(dust_example("sir"), sir)
  expect_identical(dust_example("sirode"), sirode)
  expect_identical(dust_example("walk"), walk)
  expect_error(dust_example("other"),
               "'name' must be one of 'sir', 'sirode', 'walk'")
})
