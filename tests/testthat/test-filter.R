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
    obj <- dust2_cpu_sir_alloc(base, time_start, dt, 1, 0, NULL, TRUE)
    ptr <- obj[[1]]
    dust2_cpu_sir_set_state_initial(ptr)
    incidence <- numeric(length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust2_cpu_sir_run_steps(ptr, round((time[i] - time0[i]) / dt))
      incidence[i] <- dust2_cpu_sir_state(ptr, FALSE)[5, , drop = TRUE]
    }
    sum(dpois(1:4, incidence + 1e-6, log = TRUE))
  }

  obj <- dust2_cpu_sir_unfilter_alloc(base, time_start, time, dt, data, 1, 0)
  ptr <- obj[[1]]
  expect_equal(dust2_cpu_sir_unfilter_run(ptr, NULL, NULL, FALSE), f(pars1))

  expect_equal(dust2_cpu_sir_unfilter_run(ptr, pars1, NULL, FALSE), f(pars1))
  expect_equal(dust2_cpu_sir_unfilter_run(ptr, pars2, NULL, FALSE), f(pars2))
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
    obj <- dust2_cpu_sir_alloc(pars, time_start, dt, 1, 0, NULL, TRUE)
    ptr <- obj[[1]]
    dust2_cpu_sir_set_state(ptr, state, FALSE)
    incidence <- numeric(length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust2_cpu_sir_run_steps(ptr, round((time[i] - time0[i]) / dt))
      incidence[i] <- dust2_cpu_sir_state(ptr, FALSE)[5, , drop = TRUE]
    }
    sum(dpois(1:4, incidence + 1e-6, log = TRUE))
  }

  obj <- dust2_cpu_sir_unfilter_alloc(pars, time_start, time, dt, data, 1, 0)
  ptr <- obj[[1]]
  expect_equal(dust2_cpu_sir_unfilter_run(ptr, NULL, state, FALSE), f(pars))
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
    obj <- dust2_cpu_sir_alloc(pars, time_start, dt, 1, n_groups, NULL, TRUE)
    ptr <- obj[[1]]
    dust2_cpu_sir_set_state_initial(ptr)
    incidence <- matrix(0, n_groups, length(time))
    time0 <- c(time_start, time)
    for (i in seq_along(time)) {
      dust2_cpu_sir_run_steps(ptr, round((time[i] - time0[i]) / dt))
      incidence[, i] <- dust2_cpu_sir_state(ptr, FALSE)[5, , drop = TRUE]
    }
    observed <- matrix(unlist(data, use.names = FALSE), n_groups)
    rowSums(dpois(observed, incidence + 1e-6, log = TRUE))
  }

  obj <- dust2_cpu_sir_unfilter_alloc(pars, time_start, time, dt, data, 1, 3)
  ptr <- obj[[1]]
  expect_equal(dust2_cpu_sir_unfilter_run(ptr, NULL, NULL, TRUE), f(pars))
})


test_that("validate time for filter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  time <- as.integer(c(4, 8, 12, 16))
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1

  expect_error(
    dust2_cpu_sir_unfilter_alloc(pars, 5, time, dt, data, 1, 0),
    "Expected 'time[1]' (5) to be larger than the previous value (4)",
    fixed = TRUE)
  time2 <- time + c(0, 0, .1, 0)
  expect_error(
    dust2_cpu_sir_unfilter_alloc(pars, 0, time2, dt, data, 1, 0),
    "Expected 'time[3]' to be integer-like",
    fixed = TRUE)
  expect_error(
    dust2_cpu_sir_unfilter_alloc(pars, 0, as.character(time), dt, data, 1, 0),
    "'time' must be a numeric vector",
    fixed = TRUE)
})


test_that("can run replicated unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  obj1 <- dust2_cpu_sir_unfilter_alloc(pars, time_start, time, dt, data, 5, 0)
  obj2 <- dust2_cpu_sir_unfilter_alloc(pars, time_start, time, dt, data, 1, 0)

  expect_equal(
    dust2_cpu_sir_unfilter_run(obj1[[1]], NULL, NULL, FALSE),
    rep(dust2_cpu_sir_unfilter_run(obj2[[1]], NULL, NULL, FALSE), 5))
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
  obj1 <- dust2_cpu_sir_unfilter_alloc(pars, time_start, time, dt, data, 5, 2)
  obj2 <- dust2_cpu_sir_unfilter_alloc(pars, time_start, time, dt, data, 1, 2)

  cmp <- dust2_cpu_sir_unfilter_run(obj2[[1]], NULL, NULL, TRUE)
  expect_equal(
    dust2_cpu_sir_unfilter_run(obj1[[1]], NULL, NULL, TRUE),
    matrix(rep(cmp, each = 5), 5))
})
