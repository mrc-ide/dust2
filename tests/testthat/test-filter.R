test_that("can run particle filter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), time_start, time, data,
                            n_particles = n_particles, seed = seed)
  res <- replicate(20, dust_filter_run(obj, pars))

  cmp_filter <- sir_filter_manual(
    pars, time_start, time, dt, data, n_particles, seed)
  expect_equal(res, replicate(20, cmp_filter(NULL)$log_likelihood))
})


test_that("can run particle filter and save history", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), time_start, time, data,
                            n_particles = n_particles, seed = seed)
  res1 <- dust_filter_run(obj, pars)
  expect_error(dust_filter_last_history(obj), "History is not current")
  res2 <- dust_filter_run(obj, NULL, save_history = TRUE)
  h2 <- dust_filter_last_history(obj)
  expect_equal(dim(h2), c(5, 100, 4))
  res3 <- dust_filter_run(obj, NULL)
  expect_error(dust_filter_last_history(obj), "History is not current")

  cmp_filter <- sir_filter_manual(
    pars, time_start, time, dt, data, n_particles, seed)
  cmp1 <- cmp_filter(NULL)
  cmp2 <- cmp_filter(NULL, save_history = TRUE)
  cmp3 <- cmp_filter(NULL)

  expect_equal(res1, cmp1$log_likelihood)
  expect_equal(res2, cmp2$log_likelihood)
  expect_equal(res3, cmp3$log_likelihood)
  expect_equal(h2, cmp2$history)
})


test_that("can get partial filter history", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj1 <- dust_filter_create(sir(), time_start, time, data,
                             n_particles = n_particles, seed = seed)
  obj2 <- dust_filter_create(sir(), time_start, time, data,
                             n_particles = n_particles, seed = seed,
                             index = c(2, 4))

  expect_equal(dust_filter_run(obj1, pars, save_history = TRUE),
               dust_filter_run(obj2, pars, save_history = TRUE))

  h1 <- dust_filter_last_history(obj1)
  h2 <- dust_filter_last_history(obj2)
  expect_equal(dim(h1), c(5, 100, 4))
  expect_equal(dim(h2), c(2, 100, 4))
  expect_equal(h2, h1[c(2, 4), , , drop = FALSE])
})


test_that("can run a nested particle filter and get the same result", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6))

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) {
    list(list(incidence = i), list(incidence = i + 1))
  })
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), time_start, time, data,
                            n_particles = n_particles, n_groups = 2,
                            seed = seed)

  ## Here, we can check the layout of the rng within the filter and system:
  n_streams <- (n_particles + 1) * 2
  r <- mcstate2::mcstate_rng$new(n_streams = n_streams, seed = seed)$state()
  s <- dust_filter_rng_state(obj)
  expect_equal(s, r)

  res <- replicate(20, dust_filter_run(obj, pars, save_history = TRUE))

  ## now compare:
  data1 <- lapply(data, "[[", 1)
  obj1 <- dust_filter_create(sir(), time_start, time, data1,
                             n_particles = n_particles, seed = seed)
  s1 <- dust_filter_rng_state(obj1)
  res1 <- replicate(20, dust_filter_run(obj1, pars[[1]], save_history = TRUE))
  expect_equal(res1, res[1, ])
  expect_equal(s1, s[1:3232])

  seed2 <- r[3233:3264]
  data2 <- lapply(data, "[[", 2)
  obj2 <- dust_filter_create(sir(), time_start, time, data2,
                             n_particles = n_particles, seed = seed2)
  s2 <- dust_filter_rng_state(obj2)
  res2 <- replicate(20, dust_filter_run(obj2, pars[[2]], save_history = TRUE))
  expect_equal(res2, res[2, ])
  expect_equal(s2, s[3233:6464])

  h <- dust_filter_last_history(obj)
  expect_equal(dim(h), c(5, 100, 2, 4))
  expect_equal(h[, , 1, ], dust_filter_last_history(obj1))
  expect_equal(h[, , 2, ], dust_filter_last_history(obj2))
})


test_that("can run filter and change parameters", {
  base <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  update <- list(beta = 0.15, gamma = 0.25, I0 = 15)
  pars <- modifyList(base, update)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj1 <- dust_filter_create(sir(), time_start, time, data,
                             n_particles = n_particles, seed = seed)
  obj2 <- dust_filter_create(sir(), time_start, time, data,
                             n_particles = n_particles, seed = seed)

  expect_equal(dust_filter_run(obj1, base),
               dust_filter_run(obj2, base))
  expect_equal(dust_filter_run(obj1, update),
               dust_filter_run(obj2, pars))
})


test_that("can run particle filter with manual initial state", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  state <- matrix(c(1000 - 17, 17, 0, 0, 0), ncol = 1)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), time_start, time, data,
                            n_particles = n_particles, seed = seed)
  res <- replicate(20, dust_filter_run(obj, pars, initial = state))

  cmp_filter <- sir_filter_manual(
    pars, time_start, time, dt, data, n_particles, seed)
  expect_equal(res, replicate(20, cmp_filter(NULL, state)$log_likelihood))
})


test_that("can set rng state into the filter", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6))

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) {
    list(list(incidence = i), list(incidence = i + 1))
  })
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj1 <- dust_filter_create(sir(), time_start, time, data,
                             n_particles = n_particles, n_groups = 2,
                             seed = 42)
  obj2 <- dust_filter_create(sir(), time_start, time, data,
                             n_particles = n_particles, n_groups = 2,
                             seed = 43)
  s1 <- dust_filter_rng_state(obj1)
  s2 <- dust_filter_rng_state(obj2)
  expect_false(identical(s1, s2))

  dust_filter_set_rng_state(obj2, s1)
  expect_identical(dust_filter_rng_state(obj2), s1)
  expect_identical(dust_filter_run(obj1, pars), dust_filter_run(obj2, pars))
})
