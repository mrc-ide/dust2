test_that("can compile and run with openmp", {
  skip_on_cran()
  skip_for_compilation()

  gen <- dust_compile("examples/openmp.cpp", quiet = TRUE, debug = TRUE)
  obj <- dust_system_create(gen(), list(), n_particles = 10, seed = 1)
  dust_system_set_state_initial(obj)
  s <- dust_system_state(obj)
  has_openmp <- all(s > 0)
  if (!has_openmp) {
    testthat::skip("no openmp support")
  }

  expect_equal(dust_system_state(obj), matrix(0, 1, 10))

  obj <- dust_system_create(gen(), list(), n_particles = 10, seed = 1,
                            n_threads = 5)
  dust_system_set_state_initial(obj)
  expect_equal(dust_system_state(obj), matrix(rep(0:4, each = 2), 1))

  obj <- dust_system_create(gen(), list(list(), list(), list()),
                            n_particles = 10, n_groups = 3,
                            seed = 1, n_threads = 5)
  dust_system_set_state_initial(obj)
  expect_equal(dust_system_state(obj),
               array(rep(0:4, each = 6), c(1, 10, 3)))
})
