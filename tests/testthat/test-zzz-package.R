test_that("can generate basic test package", {
  path <- withr::local_tempfile()
  create_test_package(path)
  dust_package(path, quiet = TRUE)

  pkgbuild::compile_dll(path, quiet = TRUE)
  res <- pkgload::load_all(path, quiet = TRUE)
  w <- dust_system_create(res$env$walk(), list(sd = 1), 100)
  expect_equal(dust_system_state(w), matrix(0, 1, 100))

  rm(w)
  gc()
})
