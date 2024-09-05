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
