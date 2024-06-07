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
    obj <- dust_system_create(sir(), base, time = time_start, dt = dt,
                              n_particles = 1, deterministic = TRUE)
    dust_system_set_state_initial(obj)
    incidence <- numeric(length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust_system_run_steps(obj, round((time[i] - time0[i]) / dt))
      incidence[i] <- dust_system_state(obj)[5, , drop = TRUE]
    }
    sum(dpois(1:4, incidence + 1e-6, log = TRUE))
  }

  obj <- dust_unfilter_create(sir(), time_start, time, data)
  expect_equal(dust_unfilter_run(obj, base), f(pars1))

  expect_equal(dust_unfilter_run(obj, pars = pars1), f(pars1))
  expect_equal(dust_unfilter_run(obj, pars = pars2), f(pars2))
})


test_that("can get unfilter history", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1

  obj <- dust_unfilter_create(sir(), time_start, time, data)
  dust_unfilter_run(obj, pars)
  expect_error(
    dust_unfilter_last_history(obj),
    "History is not current")
  dust_unfilter_run(obj, NULL, save_history = TRUE)
  h <- dust_unfilter_last_history(obj)
  expect_equal(dust_unfilter_last_history(obj), h)
  dust_unfilter_run(obj, NULL, save_history = FALSE)
  expect_error(
    dust_unfilter_last_history(obj),
    "History is not current")

  m <- dust_system_create(sir(), pars, time = time_start, n_particles = 1,
                          deterministic = TRUE)
  dust_system_set_state_initial(m)
  cmp <- dust_system_simulate(m, time)
  expect_equal(h, cmp)
})


test_that("can get partial unfilter history", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1

  obj1 <- dust_unfilter_create(sir(), time_start, time, data)
  obj2 <- dust_unfilter_create(sir(), time_start, time, data,
                               index = c(2, 4))
  expect_equal(dust_unfilter_run(obj1, pars, save_history = TRUE),
               dust_unfilter_run(obj2, pars, save_history = TRUE))

  h1 <- dust_unfilter_last_history(obj1)
  h2 <- dust_unfilter_last_history(obj2)
  expect_equal(dim(h1), c(5, 1, 4))
  expect_equal(dim(h2), c(2, 1, 4))
  expect_equal(h2, h1[c(2, 4), , , drop = FALSE])
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
    obj <- dust_system_create(sir(), pars, time = time_start, dt = dt,
                              n_particles = 1, deterministic = TRUE)
    dust_system_set_state(obj, state)
    incidence <- numeric(length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust_system_run_steps(obj, round((time[i] - time0[i]) / dt))
      incidence[i] <- dust_system_state(obj)[5, , drop = TRUE]
    }
    sum(dpois(1:4, incidence + 1e-6, log = TRUE))
  }

  obj <- dust_unfilter_create(sir(), time_start, time, data)
  expect_equal(dust_unfilter_run(obj, pars, initial = state), f(pars))
})


test_that("can run unfilter on structured system", {
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
    obj <- dust_system_create(sir(), pars, time = time_start, dt = dt,
                              n_particles = 1, n_groups = n_groups,
                              deterministic = TRUE)
    dust_system_set_state_initial(obj)
    incidence <- matrix(0, n_groups, length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust_system_run_steps(obj, round((time[i] - time0[i]) / dt))
      incidence[, i] <- dust_system_state(obj)[5, , , drop = TRUE]
    }
    observed <- matrix(unlist(data, use.names = FALSE), n_groups)
    rowSums(dpois(observed, incidence + 1e-6, log = TRUE))
  }

  obj <- dust_unfilter_create(sir(), time_start, time, data, n_groups = 3)
  expect_equal(dust_unfilter_run(obj, pars), f(pars))
})


test_that("validate time for unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  time <- as.integer(c(4, 8, 12, 16))
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1

  err <- expect_error(
    dust_unfilter_create(sir(), time_start = 5, time = time, data = data),
    "Time sequence is not strictly increasing",
    fixed = TRUE)
  expect_equal(
    err$body,
    c(x = "'time[1]' (4) must be greater than 'time_start' (5)"))

  err <- expect_error(
    dust_unfilter_create(sir(), time_start = 5, time = rev(time), data = data),
    "Time sequence is not strictly increasing",
    fixed = TRUE)
  expect_equal(
    err$body,
    c(x = "'time[2]' (12) must be greater than 'time[1]' (16)",
      x = "'time[3]' (8) must be greater than 'time[2]' (12)",
      x = "'time[4]' (4) must be greater than 'time[3]' (8)"))

  err <- expect_error(
    dust_unfilter_create(sir(), time_start = 5, time = 10:1, data = data),
    "Time sequence is not strictly increasing",
    fixed = TRUE)
  expect_equal(
    err$body,
    c(x = "'time[2]' (9) must be greater than 'time[1]' (10)",
      x = "'time[3]' (8) must be greater than 'time[2]' (9)",
      x = "'time[4]' (7) must be greater than 'time[3]' (8)",
      x = "'time[5]' (6) must be greater than 'time[4]' (7)",
      x = "...and 5 other errors"))

  time2 <- time + c(0, 0, .1, 0)
  expect_error(
    dust_unfilter_create(sir(), time_start = 0, time = time2,
                         data = data),
    "Expected 'time' to be integer",
    fixed = TRUE)
  expect_error(
    dust_unfilter_create(sir(), time_start = 0, time = as.character(time),
                         data = data),
    "Expected 'time' to be integer",
    fixed = TRUE)
})


test_that("can run replicated unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  obj1 <- dust_unfilter_create(sir(), time_start, time, data, n_particles = 5)
  obj2 <- dust_unfilter_create(sir(), time_start, time, data, n_particles = 1)
  expect_equal(
    dust_unfilter_run(obj1, pars),
    rep(dust_unfilter_run(obj2, pars), 5))
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
  obj1 <- dust_unfilter_create(sir(), time_start, time, data,
                               n_particles = 5, n_groups = 2)
  obj2 <- dust_unfilter_create(sir(), time_start, time, data,
                               n_particles = 1, n_groups = 2)
  expect_equal(
    dust_unfilter_run(obj1, pars),
    matrix(rep(dust_unfilter_run(obj2, pars), each = 5), 5))
})


test_that("can save history from structured unfilter", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6))
  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) {
    lapply(seq_len(2), function(j) list(incidence = 2 * (i - 1) + j))
  })
  data1 <- lapply(data, function(x) x[[1]])
  data2 <- lapply(data, function(x) x[[2]])
  dt <- 1
  obj <- dust_unfilter_create(sir(), time_start, time, data,
                              n_groups = 2)
  obj1 <- dust_unfilter_create(sir(), time_start, time, data1)
  obj2 <- dust_unfilter_create(sir(), time_start, time, data2)

  ll <- dust_unfilter_run(obj, pars, save_history = TRUE)
  ll1 <- dust_unfilter_run(obj1, pars[[1]], save_history = TRUE)
  ll2 <- dust_unfilter_run(obj2, pars[[2]], save_history = TRUE)

  expect_equal(ll, c(ll1, ll2))
  h <- dust_unfilter_last_history(obj)
  h1 <- dust_unfilter_last_history(obj1)
  h2 <- dust_unfilter_last_history(obj2)

  expect_equal(dim(h), c(5, 1, 2, 4))
  expect_equal(dim(h1), c(5, 1, 4))
  expect_equal(dim(h2), c(5, 1, 4))
  expect_equal(array(h[, , 1, ], dim(h1)), h1)
  expect_equal(array(h[, , 2, ], dim(h2)), h2)
})
