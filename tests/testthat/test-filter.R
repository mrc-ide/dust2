test_that("can run an unfilter", {
  base <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  pars1 <- list(beta = 0.1, gamma = 0.2, I0 = 10)
  pars2 <- list(beta = 0.2, gamma = 0.2, I0 = 10)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1

  ## Manually compute likelihood:
  f <- function(pars) {
    base[names(pars)] <- pars
    obj <- dust_model_create(sir(), base, time = time_start, dt = dt,
                             n_particles = 1, deterministic = TRUE)
    dust_model_set_state_initial(obj)
    incidence <- numeric(length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust_model_run_steps(obj, round((time[i] - time0[i]) / dt))
      incidence[i] <- dust_model_state(obj)[5, , drop = TRUE]
    }
    sum(dpois(1:4, incidence + 1e-6, log = TRUE))
  }

  obj <- dust_unfilter_create(sir(), base, time_start, time, data)
  expect_equal(dust_unfilter_run(obj), f(pars1))

  expect_equal(dust_unfilter_run(obj, pars = pars1), f(pars1))
  expect_equal(dust_unfilter_run(obj, pars = pars2), f(pars2))
})


test_that("can run an unfilter with manually set state", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  state <- matrix(c(1000 - 17, 17, 0, 0, 0), ncol = 1)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1

  ## Manually compute likelihood:
  f <- function(pars) {
    obj <- dust_model_create(sir(), pars, time = time_start, dt = dt,
                             n_particles = 1, deterministic = TRUE)
    dust_model_set_state(obj, state)
    incidence <- numeric(length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust_model_run_steps(obj, round((time[i] - time0[i]) / dt))
      incidence[i] <- dust_model_state(obj)[5, , drop = TRUE]
    }
    sum(dpois(1:4, incidence + 1e-6, log = TRUE))
  }

  obj <- dust_unfilter_create(sir(), pars, time_start, time, data)
  expect_equal(dust_unfilter_run(obj, initial = state), f(pars))
})


test_that("can run unfilter on structured model", {
  base <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  n_groups <- 3
  pars <- lapply(seq_len(n_groups),
                 function(i) modifyList(base, list(beta = i * 0.1)))

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) {
    lapply(seq_len(n_groups), function(j) list(incidence = 2 * (i - 1) + j))
  })
  dt <- 1

  ## Manually compute likelihood:
  f <- function(pars) {
    obj <- dust_model_create(sir(), pars, time = time_start, dt = dt,
                             n_particles = 1, n_groups = n_groups,
                             deterministic = TRUE)
    dust_model_set_state_initial(obj)
    incidence <- matrix(0, n_groups, length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust_model_run_steps(obj, round((time[i] - time0[i]) / dt))
      incidence[, i] <- dust_model_state(obj)[5, , , drop = TRUE]
    }
    observed <- matrix(unlist(data, use.names = FALSE), n_groups)
    rowSums(dpois(observed, incidence + 1e-6, log = TRUE))
  }

  obj <- dust_unfilter_create(sir(), pars, time_start, time, data, n_groups = 3)
  expect_equal(dust_unfilter_run(obj), f(pars))
})


test_that("validate time for filter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  time <- as.integer(c(4, 8, 12, 16))
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1

  expect_error(
    dust_unfilter_create(sir(), pars, time_start = 5, time = time, data = data),
    "Expected 'time[1]' (5) to be larger than the previous value (4)",
    fixed = TRUE)
  time2 <- time + c(0, 0, .1, 0)
  expect_error(
    dust_unfilter_create(sir(), pars, time_start = 0, time = time2,
                         data = data),
    "Expected 'time[3]' to be integer-like",
    fixed = TRUE)
  expect_error(
    dust_unfilter_create(sir(), pars, time_start = 0, time = as.character(time),
                         data = data),
    "'time' must be a numeric vector",
    fixed = TRUE)
})


test_that("can run replicated unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  obj1 <- dust_unfilter_create(sir(), pars, time_start, time, data,
                               n_particles = 5)
  obj2 <- dust_unfilter_create(sir(), pars, time_start, time, data,
                               n_particles = 1)
  expect_equal(
    dust_unfilter_run(obj1),
    rep(dust_unfilter_run(obj2), 5))
})


test_that("can run replicated structured unfilter", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6))
  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) {
    lapply(seq_len(2), function(j) list(incidence = 2 * (i - 1) + j))
  })
  dt <- 1
  obj1 <- dust_unfilter_create(sir(), pars, time_start, time, data,
                               n_particles = 5, n_groups = 2)
  obj2 <- dust_unfilter_create(sir(), pars, time_start, time, data,
                               n_particles = 1, n_groups = 2)
  expect_equal(
    dust_unfilter_run(obj1),
    matrix(rep(dust_unfilter_run(obj2), each = 5), 5))
})


test_that("can run particle filter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), pars, time_start, time, data,
                           n_particles = n_particles, seed = seed)
  res <- replicate(20, dust_filter_run(obj))

  cmp_filter <- sir_filter_manual(
    pars, time_start, time, dt, data, n_particles, seed)
  expect_equal(res, replicate(20, cmp_filter(NULL)))
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

  obj <- dust_filter_create(sir(), pars, time_start, time, data,
                            n_particles = n_particles, n_groups = 2,
                            seed = seed)

  ## Here, we can check the layout of the rng within the filter and model:
  n_streams <- (n_particles + 1) * 2
  r <- mcstate2::mcstate_rng$new(n_streams = n_streams, seed = seed)$state()
  s <- dust_filter_rng_state(obj)
  expect_equal(s, r)

  res <- replicate(20, dust_filter_run(obj))

  ## now compare:
  data1 <- lapply(data, "[[", 1)
  obj1 <- dust_filter_create(sir(), pars[[1]], time_start, time, data1,
                             n_particles = n_particles, seed = seed)
  s1 <- dust_filter_rng_state(obj1)
  res1 <- replicate(20, dust_filter_run(obj1))
  expect_equal(res1, res[1, ])
  expect_equal(s1, s[1:3232])

  seed2 <- r[3233:3264]
  data2 <- lapply(data, "[[", 2)
  obj2 <- dust_filter_create(sir(), pars[[2]], time_start, time, data2,
                             n_particles = n_particles, seed = seed2)
  s2 <- dust_filter_rng_state(obj2)
  res2 <- replicate(20, dust_filter_run(obj2))
  expect_equal(res2, res[2, ])
  expect_equal(s2, s[3233:6464])
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

  obj1 <- dust_filter_create(sir(), base, time_start, time, data,
                             n_particles = n_particles, seed = seed)
  obj2 <- dust_filter_create(sir(), pars, time_start, time, data,
                             n_particles = n_particles, seed = seed)

  expect_equal(dust_filter_run(obj1, update),
               dust_filter_run(obj2))
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

  obj <- dust_filter_create(sir(), pars, time_start, time, data,
                           n_particles = n_particles, seed = seed)
  res <- replicate(20, dust_filter_run(obj, initial = state))

  cmp_filter <- sir_filter_manual(
    pars, time_start, time, dt, data, n_particles, seed)
  expect_equal(res, replicate(20, cmp_filter(NULL, state)))
})
