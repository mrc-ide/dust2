test_that("can run an unfilter", {
  base <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  pars1 <- list(beta = 0.1, gamma = 0.2, I0 = 10)
  pars2 <- list(beta = 0.2, gamma = 0.2, I0 = 10)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  dt <- 1

  ## Manually compute likelihood:
  f <- function(pars) {
    base[names(pars)] <- pars
    obj <- dust_system_create(sir(), base, time = time_start, dt = dt,
                              n_particles = 1, deterministic = TRUE)
    dust_system_set_state_initial(obj)
    incidence <- numeric(nrow(data))
    time0 <- c(time_start, data$time)
    for (i in seq_len(nrow(data))) {
      dust_system_run_to_time(obj, data$time[i])
      incidence[i] <- dust_system_state(obj)[5]
    }
    sum(dpois(data$incidence, incidence + 1e-6, log = TRUE))
  }

  obj <- dust_unfilter_create(sir(), time_start, data)
  expect_equal(dust_unfilter_run(obj, base), f(pars1))

  expect_equal(dust_unfilter_run(obj, pars = pars1), f(pars1))
  expect_equal(dust_unfilter_run(obj, pars = pars2), f(pars2))
})


test_that("can only use pars = NULL on initialised unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sir(), time_start, data)
  expect_error(dust_unfilter_run(obj, NULL),
               "'pars' cannot be NULL, as unfilter is not initialised")
  ll <- dust_unfilter_run(obj, pars)
  expect_identical(dust_unfilter_run(obj, NULL), ll)
})


test_that("can get unfilter history", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sir(), time_start, data)
  expect_error(
    dust_unfilter_last_history(obj),
    "History is not current")
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
  cmp <- dust_system_simulate(m, data$time)
  expect_equal(h, cmp)
})


test_that("can get partial unfilter history", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj1 <- dust_unfilter_create(sir(), time_start, data)
  obj2 <- dust_unfilter_create(sir(), time_start, data,
                               index_state = c(2, 4))
  expect_equal(dust_unfilter_run(obj1, pars, save_history = TRUE),
               dust_unfilter_run(obj2, pars, save_history = TRUE))

  h1 <- dust_unfilter_last_history(obj1)
  h2 <- dust_unfilter_last_history(obj2)
  expect_equal(dim(h1), c(5, 4))
  expect_equal(dim(h2), c(2, 4))
  expect_equal(h2, h1[c(2, 4), , drop = FALSE])
})


test_that("can run an unfilter with manually set state", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  state <- matrix(c(1000 - 17, 17, 0, 0, 0), ncol = 1)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  dt <- 1

  ## Manually compute likelihood:
  f <- function(pars) {
    obj <- dust_system_create(sir(), pars, time = time_start, dt = dt,
                              n_particles = 1, deterministic = TRUE)
    dust_system_set_state(obj, state)
    incidence <- numeric(nrow(data))
    time0 <- c(time_start, data$time)
    for (i in seq_along(incidence)) {
      dust_system_run_to_time(obj, data$time[i])
      incidence[i] <- dust_system_state(obj)[5]
    }
    sum(dpois(data$incidence, incidence + 1e-6, log = TRUE))
  }

  obj <- dust_unfilter_create(sir(), time_start, data)
  expect_equal(dust_unfilter_run(obj, pars, initial = state), f(pars))
})


test_that("can run unfilter on structured system", {
  base <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  n_groups <- 3
  pars <- lapply(seq_len(n_groups),
                 function(i) modifyList(base, list(beta = i * 0.1)))

  time_start <- 0
  data <- data.frame(time = rep(c(4, 8, 12, 16), 3),
                     group = rep(1:3, each = 4),
                     incidence = 1:12)

  obj <- dust_unfilter_create(sir(), time_start, data, n_groups = 3)
  ll <- dust_unfilter_run(obj, pars)

  obj1 <- dust_unfilter_create(sir(), time_start, data[data$group == 1, ])
  obj2 <- dust_unfilter_create(sir(), time_start, data[data$group == 2, ])
  obj3 <- dust_unfilter_create(sir(), time_start, data[data$group == 3, ])
  ll1 <- dust_unfilter_run(obj1, pars[[1]])
  ll2 <- dust_unfilter_run(obj2, pars[[2]])
  ll3 <- dust_unfilter_run(obj3, pars[[3]])

  expect_equal(ll, c(ll1, ll2, ll3))
})


test_that("can run replicated unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj1 <- dust_unfilter_create(sir(), time_start, data, n_particles = 5)
  obj2 <- dust_unfilter_create(sir(), time_start, data, n_particles = 1)
  expect_equal(
    dust_unfilter_run(obj1, pars),
    rep(dust_unfilter_run(obj2, pars), 5))
})


test_that("can run replicated structured unfilter", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6))
  time_start <- 0
  data <- data.frame(time = rep(c(4, 8, 12, 16), 2),
                     group = rep(1:2, each = 4),
                     incidence = 1:8)
  obj1 <- dust_unfilter_create(sir(), time_start, data,
                               n_particles = 5, n_groups = 2)
  obj2 <- dust_unfilter_create(sir(), time_start, data,
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
  data <- data.frame(time = rep(c(4, 8, 12, 16), 2),
                     group = rep(1:2, each = 4),
                     incidence = 1:8)
  data1 <- data[data$group == 1, -2]
  data2 <- data[data$group == 2, -2]
  obj <- dust_unfilter_create(sir(), time_start, data)
  obj1 <- dust_unfilter_create(sir(), time_start, data1)
  obj2 <- dust_unfilter_create(sir(), time_start, data2)

  ll <- dust_unfilter_run(obj, pars, save_history = TRUE)
  ll1 <- dust_unfilter_run(obj1, pars[[1]], save_history = TRUE)
  ll2 <- dust_unfilter_run(obj2, pars[[2]], save_history = TRUE)

  expect_equal(ll, c(ll1, ll2))
  h <- dust_unfilter_last_history(obj)
  h1 <- dust_unfilter_last_history(obj1)
  h2 <- dust_unfilter_last_history(obj2)

  expect_equal(dim(h), c(5, 2, 4))
  expect_equal(dim(h1), c(5, 4))
  expect_equal(dim(h2), c(5, 4))
  expect_equal(array(h[, 1, ], dim(h1)), h1)
  expect_equal(array(h[, 2, ], dim(h2)), h2)
})


test_that("history output can have dimensions preserved", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj1 <- dust_unfilter_create(sir(), time_start, data)
  ll1 <- dust_unfilter_run(obj1, pars, save_history = TRUE)

  obj2 <- dust_unfilter_create(sir(), time_start, data,
                               preserve_particle_dimension = TRUE)
  ll2 <- dust_unfilter_run(obj2, pars, save_history = TRUE)
  expect_equal(ll2, ll1)

  obj3 <- dust_unfilter_create(sir(), time_start, data,
                               preserve_group_dimension = TRUE)
  ll3 <- dust_unfilter_run(obj3, list(pars), save_history = TRUE)
  expect_equal(ll3, ll1)

  obj4 <- dust_unfilter_create(sir(), time_start, data,
                               preserve_group_dimension = TRUE,
                               preserve_particle_dimension = TRUE)
  ll4 <- dust_unfilter_run(obj4, list(pars), save_history = TRUE)
  expect_equal(ll4, matrix(ll1, 1, 1))

  h1 <- dust_unfilter_last_history(obj1)
  expect_equal(dim(h1), c(5, 4))

  h2 <- dust_unfilter_last_history(obj2)
  expect_equal(dim(h2), c(5, 1, 4))
  expect_equal(drop(h2), h1)

  h3 <- dust_unfilter_last_history(obj3)
  expect_equal(dim(h3), c(5, 1, 4))
  expect_equal(drop(h3), h1)

  h4 <- dust_unfilter_last_history(obj4)
  expect_equal(dim(h4), c(5, 1, 1, 4))
  expect_equal(drop(h4), h1)
})


test_that("can extract a subset of an unfilter run in its entirety", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6))
  time_start <- 0
  data <- data.frame(time = rep(c(4, 8, 12, 16), 2),
                     group = rep(1:2, each = 4),
                     incidence = 1:8)
  obj <- dust_unfilter_create(sir(), time_start, data, n_groups = 2)
  dust_unfilter_run(obj, pars, save_history = TRUE)
  h1 <- dust_unfilter_last_history(obj, index_group = 1)
  h2 <- dust_unfilter_last_history(obj, index_group = 2)
  h <- dust_unfilter_last_history(obj)
  expect_equal(h1, h[, 1, , drop = FALSE])
  expect_equal(h2, h[, 2, , drop = FALSE])
})


test_that("can extract a state-filtered subset of an unfilter run", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6))
  time_start <- 0
  data <- data.frame(time = rep(c(4, 8, 12, 16), 2),
                     group = rep(1:2, each = 4),
                     incidence = 1:8)
  obj <- dust_unfilter_create(sir(), time_start, data, n_groups = 2,
                              index_state = c(1, 3, 5))
  dust_unfilter_run(obj, pars, save_history = TRUE)
  h1 <- dust_unfilter_last_history(obj, index_group = 1)
  h2 <- dust_unfilter_last_history(obj, index_group = 2)
  h <- dust_unfilter_last_history(obj)
  expect_equal(h1, h[, 1, , drop = FALSE])
  expect_equal(h2, h[, 2, , drop = FALSE])
})


test_that("can run a subset of an unfilter", {
  pars0 <- rep(
    list(list(beta = 0.01, gamma = 0.01, N = 100, I0 = 1, exp_noise = 1e6)),
    3)
  pars1 <- list(
    list(beta = 0.1, gamma = 0.2, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, I0 = 10, exp_noise = 1e6),
    list(beta = 0.3, gamma = 0.2, I0 = 10, exp_noise = 1e6))

  time_start <- 0
  data <- data.frame(time = rep(c(4, 8, 12, 16), 3),
                     group = rep(1:3, each = 4),
                     incidence = 1:12)

  obj1 <- dust_unfilter_create(sir(), time_start, data)
  ll0 <- dust_unfilter_run(obj1, pars0)
  ll1 <- dust_unfilter_run(obj1, pars1, save_history = TRUE)
  expect_length(ll1, 3)
  h1 <- dust_unfilter_last_history(obj1)

  obj2 <- dust_unfilter_create(sir(), time_start, data)
  ll0 <- dust_unfilter_run(obj2, pars0)
  expect_equal(
    dust_unfilter_run(obj2, pars1[3:1], index_group = 3:1, save_history = TRUE),
    rev(ll1))

  expect_equal(dust_unfilter_last_history(obj2, index_group = 3:1),
               h1[, 3:1, ])

  for (i in 1:3) {
    ll_i <- dust_unfilter_run(obj2, pars1[i], index_group = i,
                              save_history = TRUE)
    expect_equal(dust_unfilter_last_history(obj2, index_group = i),
                 h1[, i, , drop = FALSE])
    j <- i %% 3 + 1
    expect_error(dust_unfilter_last_history(obj2, index_group = j),
                 sprintf("History for group '%d' is not current", j))
  }
})


test_that("can print an unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_unfilter_create(sir(), time_start, data, n_particles = 5)
  res <- evaluate_promise(withVisible(print(obj)))
  expect_mapequal(res$result,
                  list(value = obj, visible = FALSE))
  expect_match(res$messages, "<dust_unfilter (sir)>",
               fixed = TRUE, all = FALSE)
  expect_match(res$messages, "5 particles", all = FALSE)
})


test_that("can extract final state from an unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_unfilter_create(sir(), time_start, data)

  ## Looks like I have broken the _old_ history saving somehow?
  dust_unfilter_run(obj, pars, save_history = TRUE)
  h <- dust_unfilter_last_history(obj)
  s <- dust_unfilter_last_state(obj)
  expect_equal(s, h[, 4])

  dust_unfilter_run(obj, pars, save_history = FALSE)
  expect_error(dust_unfilter_last_history(obj), "History is not current")
  expect_equal(dust_unfilter_last_state(obj), s)
})


test_that("can extract final state from an unfilter, ignoring state index", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_unfilter_create(sir(), time_start, data, index_state = 1:3)

  dust_unfilter_run(obj, pars, save_history = TRUE)
  h <- dust_unfilter_last_history(obj)
  s <- dust_unfilter_last_state(obj)
  expect_equal(s[1:3], h[, 4])
  expect_length(s, 5)
})


test_that("can extract final state from grouped unfilter", {
  pars <- list(
    list(beta = 0.1, gamma = 0.2, I0 = 10, exp_noise = 1e6),
    list(beta = 0.2, gamma = 0.2, I0 = 10, exp_noise = 1e6),
    list(beta = 0.3, gamma = 0.2, I0 = 10, exp_noise = 1e6))

  time_start <- 0
  data <- data.frame(time = rep(c(4, 8, 12, 16), 3),
                     group = rep(1:3, each = 4),
                     incidence = 1:12)

  obj <- dust_unfilter_create(sir(), time_start, data)
  dust_unfilter_run(obj, pars, save_history = TRUE)
  h <- dust_unfilter_last_history(obj)
  s <- dust_unfilter_last_state(obj)
  expect_equal(s, h[, , 4])

  dust_unfilter_run(obj, pars, save_history = FALSE)
  expect_error(dust_unfilter_last_history(obj), "History is not current")
  expect_equal(dust_unfilter_last_state(obj), s)
})
