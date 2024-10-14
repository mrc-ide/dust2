test_that("can use two zero_every", {
  skip_for_compilation()
  gen <- dust_compile("examples/zerotwice.cpp", quiet = TRUE, debug = TRUE)

  sys <- dust_system_create(gen(), list(a = 2, b = 300), 1)
  y <- dust_system_simulate(sys, 1:20)
  expect_equal(y[1, ], rep(1:2, length.out = 20))
  expect_equal(y[2, ], 1:20)

  sys <- dust_system_create(gen(), list(a = 300, b = 5), 1)
  y <- dust_system_simulate(sys, 1:20)
  expect_equal(y[1, ], 1:20)
  expect_equal(y[2, ], rep(1:5, length.out = 20))

  sys <- dust_system_create(gen(), list(a = 2, b = 2), 1)
  y <- dust_system_simulate(sys, 1:20)
  expect_equal(y[1, ], rep(1:2, length.out = 20))
  expect_equal(y[2, ], y[1, ])

  sys <- dust_system_create(gen(), list(a = 2, b = 5), 1)
  y <- dust_system_simulate(sys, 1:20)
  expect_equal(y[1, ], rep(1:2, length.out = 20))
  expect_equal(y[2, ], rep(1:5, length.out = 20))
})


test_that("can run a model with internal storage", {
  skip_for_compilation()
  gen <- dust_compile("examples/internal.cpp", quiet = TRUE, debug = TRUE)

  len <- 5
  set.seed(1)
  sys <- dust_system_create(gen(), list(sd = 1, len = len), 1)
  seed <- dust_system_rng_state(sys)
  y <- drop(dust_system_simulate(sys, 1:20))

  r <- monty::monty_rng$new(seed = seed)$normal(20 * len, 0, 1)
  expect_equal(y, cumsum(colMeans(matrix(r, len, 20))))
})
