test_that("Resampling works as expected", {
  ref_resample_weight <- function(w, u) {
    n <- length(w)
    uu <- u / n + seq(0, by = 1 / n, length.out = n)
    cw <- cumsum(w / sum(w))
    findInterval(uu, cw) + 1L
  }

  set.seed(1)
  w <- abs(rcauchy(20))
  u <- runif(1)

  expect_equal(
    test_resample_weight(w, u) + 1L,
    ref_resample_weight(w, u))
  expect_equal(
    test_resample_weight(w, 0) + 1L,
    ref_resample_weight(w, 0))
  expect_equal(
    test_resample_weight(w, 1 - 1e-8) + 1L,
    ref_resample_weight(w, 1 - 1e-8))
})


test_that("can use history", {
  time <- seq(0, 10, length.out = 11)
  n_state <- 6
  n_particles <- 7
  n_groups <- 3
  n_time <- length(time)
  s <- lapply(seq_along(time), function(i) {
    array(runif(n_state * n_particles * n_groups),
          c(n_state, n_particles, n_groups))
  })
  s_arr <- array(unlist(s), c(n_state, n_particles, n_groups, n_time))
  expect_equal(test_history(time, s, reorder = TRUE),
               list(time, s_arr))
  expect_equal(test_history(time, s, reorder = FALSE),
               list(time, s_arr))
  expect_equal(test_history(time, s[1:3], reorder = TRUE),
               list(time[1:3], s_arr[, , , 1:3]))
  expect_equal(test_history(time, s, order = vector("list", length(time)),
                            reorder = TRUE),
               list(time, s_arr))
})


test_that("can scale log weights", {
  w <- log(abs(rcauchy(20)))
  expect_equal(test_scale_log_weights(w[1]), list(w[1], w[1]))
  expect_equal(test_scale_log_weights(w),
               list(log(mean(exp(w))), exp(w - max(w))))
  expect_equal(test_scale_log_weights(c(w, NA)),
               list(log(mean(c(exp(w), 0))), c(exp(w - max(w)), 0)))
  expect_equal(test_scale_log_weights(c(-Inf, -Inf)),
               list(-Inf, c(-Inf, -Inf)))
})


test_that("can reorder history with no groups", {
  ## This is really hard to get right so let's actually simulate forward:
  time <- seq(0, 10, length.out = 11)
  n_time <- length(time)
  n_state <- 6
  n_particles <- 7
  n_groups <- 1
  state <- vector("list", length(time))
  order <- vector("list", length(time))
  true <- array(NA_real_, c(n_state, n_particles, n_groups, n_time))
  s <- array(0, c(n_state, n_particles, n_groups))
  set.seed(1)
  for (i in seq_along(time)) {
    s <- s + runif(length(s))
    if (i > 1 && i %% 2 == 0) {
      k <- replicate(n_groups, sample(n_particles, replace = TRUE))
      order[[i]] <- as.integer(k - 1L)
      for (j in seq_len(n_groups)) {
        s[, , j] <- s[, k[, j], j] + 1
        true[, , j, seq_len(i - 1)] <- true[, k[, j], j, seq_len(i - 1)]
      }
    }
    state[[i]] <- s
    true[, , , i] <- s
  }
  state_arr <- array(unlist(state), dim(true))

  ## Pass in, but ignore index
  expect_equal(
    test_history(time, state, order = order, reorder = FALSE),
    list(time, state_arr))

  ## Really simple, add an index that does not reorder anything:
  expect_equal(
    test_history(time, state[1], order = order[1], reorder = TRUE),
    list(time[1], state_arr[, , , 1, drop = FALSE]))
  expect_equal(
    test_history(time, state[1:2], order = list(NULL, 0:6), reorder = TRUE),
    list(time[1:2], state_arr[, , , 1:2, drop = FALSE]))

  ## Proper reordering with the full index:
  expect_equal(
    test_history(time, state, order = order, reorder = TRUE),
    list(time, true))
})


test_that("can reorder history on the way out", {
  ## This is really hard to get right so let's actually simulate forward:
  time <- seq(0, 10, length.out = 11)
  n_time <- length(time)
  n_state <- 6
  n_particles <- 7
  n_groups <- 3
  state <- vector("list", length(time))
  order <- vector("list", length(time))
  true <- array(NA_real_, c(n_state, n_particles, n_groups, n_time))
  s <- array(0, c(n_state, n_particles, n_groups))
  set.seed(1)
  for (i in seq_along(time)) {
    s <- s + runif(length(s))
    if (i > 1 && i %% 2 == 0) {
      k <- replicate(n_groups, sample(n_particles, replace = TRUE))
      order[[i]] <- as.integer(k - 1L)
      for (j in seq_len(n_groups)) {
        s[, , j] <- s[, k[, j], j] + 1
        true[, , j, seq_len(i - 1)] <- true[, k[, j], j, seq_len(i - 1)]
      }
    }
    state[[i]] <- s
    true[, , , i] <- s
  }

  state_arr <- array(unlist(state), dim(true))
  expect_equal(test_history(time, state, order = order, reorder = FALSE),
               list(time, state_arr))
  expect_equal(test_history(time, state, order = order, reorder = TRUE),
               list(time, true))
})


test_that("can extract history with group index, no reordering", {
  time <- seq(0, 10, length.out = 11)
  n_state <- 5
  n_particles <- 7
  n_groups <- 3
  n_time <- length(time)
  s <- lapply(seq_along(time), function(i) {
    array(runif(n_state * n_particles * n_groups),
          c(n_state, n_particles, n_groups))
  })

  s_arr <- array(unlist(s), c(n_state, n_particles, n_groups, n_time))

  expect_equal(test_history(time, s, index_group = NULL),
               list(time, s_arr))
  expect_equal(test_history(time, s, index_group = seq_len(n_groups)),
               list(time, s_arr))

  expect_equal(test_history(time, s, index_group = 2),
               list(time, s_arr[, , 2, , drop = FALSE]))
  expect_equal(test_history(time, s, index_group = c(3, 1)),
               list(time, s_arr[, , c(3, 1), , drop = FALSE]))
})
