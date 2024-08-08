test_that("can run particle filter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), time_start, data, dt = dt,
                            n_particles = n_particles, seed = seed)
  res <- replicate(20, dust_filter_run(obj, pars))

  cmp_filter <- sir_filter_manual(
    pars, time_start, data, dt, n_particles, seed)
  expect_equal(res, replicate(20, cmp_filter(NULL)$log_likelihood))
})


test_that("can only use pars = NULL on initialised filter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj1 <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                             seed = 42)
  obj2 <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                             seed = 42)
  expect_error(dust_filter_run(obj1, NULL),
               "'pars' cannot be NULL, as filter is not initialised")
  expect_identical(dust_filter_run(obj1, pars),
                   dust_filter_run(obj2, pars))
  expect_identical(dust_filter_run(obj2, NULL),
                   dust_filter_run(obj1, pars))
})


test_that("can run particle filter and save history", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), time_start, data,
                            n_particles = n_particles, seed = seed)
  expect_error(dust_filter_last_history(obj), "History is not current")
  res1 <- dust_filter_run(obj, pars)
  expect_error(dust_filter_last_history(obj), "History is not current")
  res2 <- dust_filter_run(obj, NULL, save_history = TRUE)
  h2 <- dust_filter_last_history(obj)
  expect_equal(dim(h2), c(5, 100, 4))
  res3 <- dust_filter_run(obj, NULL)
  expect_error(dust_filter_last_history(obj), "History is not current")

  cmp_filter <- sir_filter_manual(
    pars, time_start, data, dt, n_particles, seed)
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
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  n_particles <- 100
  seed <- 42

  obj1 <- dust_filter_create(sir(), time_start, data,
                             n_particles = n_particles, seed = seed)
  obj2 <- dust_filter_create(sir(), time_start, data,
                             n_particles = n_particles, seed = seed,
                             index_state = c(2, 4))

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
  data <- data.frame(time = rep(c(4, 8, 12, 16), 2),
                     group = rep(1:2, each = 4),
                     incidence = c(1:4, 2:5))
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), time_start, data,
                            n_particles = n_particles, n_groups = 2,
                            seed = seed)

  ## Here, we can check the layout of the rng within the filter and system:
  n_streams <- (n_particles + 1) * 2
  r <- mcstate2::mcstate_rng$new(n_streams = n_streams, seed = seed)$state()
  s <- dust_filter_rng_state(obj)
  expect_equal(s, r)

  expect_true(obj$preserve_group_dimension)

  res <- replicate(20, dust_filter_run(obj, pars, save_history = TRUE))

  ## now compare
  data1 <- data[data$group == 1, -2]
  obj1 <- dust_filter_create(sir(), time_start, data1,
                             n_particles = n_particles, seed = seed)
  s1 <- dust_filter_rng_state(obj1)
  res1 <- replicate(20, dust_filter_run(obj1, pars[[1]], save_history = TRUE))
  expect_equal(res1, res[1, ])
  expect_equal(s1, s[1:3232])

  seed2 <- r[3233:3264]
  data2 <- data[data$group == 2, -2]
  obj2 <- dust_filter_create(sir(), time_start, data2,
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
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  n_particles <- 100
  seed <- 42

  obj1 <- dust_filter_create(sir(), time_start, data,
                             n_particles = n_particles, seed = seed)
  obj2 <- dust_filter_create(sir(), time_start, data,
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
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), time_start, data,
                            n_particles = n_particles, seed = seed)
  res <- replicate(20, dust_filter_run(obj, pars, initial = state))

  cmp_filter <- sir_filter_manual(
    pars, time_start, data, dt, n_particles, seed)
  expect_equal(res, replicate(20, cmp_filter(NULL, state)$log_likelihood))
})


test_that("can set rng state into the filter before running", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6))

  time_start <- 0
  data <- data.frame(time = rep(c(4, 8, 12, 16), 2),
                     group = rep(1:2, each = 4),
                     incidence = c(1:4, 2:5))
  n_particles <- 100
  seed <- 42

  obj1 <- dust_filter_create(sir(), time_start, data,
                             n_particles = n_particles, n_groups = 2,
                             seed = 42)
  obj2 <- dust_filter_create(sir(), time_start, data,
                             n_particles = n_particles, n_groups = 2,
                             seed = 43)
  s1 <- dust_filter_rng_state(obj1)
  s2 <- dust_filter_rng_state(obj2)
  expect_false(identical(s1, s2))

  dust_filter_set_rng_state(obj2, s1)
  expect_identical(dust_filter_rng_state(obj2), s1)
  expect_identical(dust_filter_run(obj1, pars), dust_filter_run(obj2, pars))
})


test_that("can set rng state into the filter after running", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6))

  time_start <- 0
  data <- data.frame(time = rep(c(4, 8, 12, 16), 2),
                     group = rep(1:2, each = 4),
                     incidence = c(1:4, 2:5))
  n_particles <- 100
  seed <- 42

  obj1 <- dust_filter_create(sir(), time_start, data,
                             n_particles = n_particles, n_groups = 2,
                             seed = 42)
  obj2 <- dust_filter_create(sir(), time_start, data,
                             n_particles = n_particles, n_groups = 2,
                             seed = 43)
  expect_false(identical(dust_filter_run(obj1, pars),
                         dust_filter_run(obj2, pars)))

  s1 <- dust_filter_rng_state(obj1)
  s2 <- dust_filter_rng_state(obj2)
  expect_false(identical(s1, s2))

  dust_filter_set_rng_state(obj2, s1)
  expect_identical(dust_filter_rng_state(obj2), s1)
  expect_identical(dust_filter_run(obj1, pars), dust_filter_run(obj2, pars))
})


test_that("can create a new filter via copying", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  n_particles <- 100
  seed <- 42

  obj1 <- dust_filter_create(sir(), time_start, data,
                             n_particles = n_particles, seed = seed)
  obj2 <- dust_filter_copy(obj1)
  obj3 <- dust_filter_copy(obj1, seed = seed)

  expect_setequal(names(obj1), names(obj2))
  nms <- setdiff(names(obj1), "initial_rng_state")
  expect_equal(as.list.environment(obj1)[nms],
               as.list.environment(obj2)[nms])
  expect_false(identical(obj1$initial_rng_state, obj2$initial_rng_state))

  expect_equal(obj1, obj3)
  expect_true(identical(obj1$initial_rng_state, obj3$initial_rng_state))

  ll1 <- dust_filter_run(obj1, pars)
  ll2 <- dust_filter_run(obj2, pars)
  ll3 <- dust_filter_run(obj3, pars)

  expect_false(identical(ll1, ll2))
  expect_identical(ll1, ll3)
})
