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


test_that("can use trajectories", {
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

  expect_equal(test_trajectories(time, s, reorder = TRUE),
               list(time, s_arr, NULL))
  expect_equal(test_trajectories(time, s, reorder = FALSE),
               list(time, s_arr, NULL))
  expect_equal(test_trajectories(time, s[1:3], reorder = TRUE),
               list(time[1:3], s_arr[, , , 1:3], NULL))
  expect_equal(test_trajectories(time, s, order = vector("list", length(time)),
                                 reorder = TRUE),
               list(time, s_arr, NULL))

  res <- test_trajectories(time, s, select_particle = c(6, 4, 2))[[2]]
  expect_equal(dim(res), c(n_state, n_groups, n_time))
  expect_equal(res[, 1, ], s_arr[, 6, 1, ])
  expect_equal(res[, 2, ], s_arr[, 4, 2, ])
  expect_equal(res[, 3, ], s_arr[, 2, 3, ])
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


test_that("can reorder trajectories with no groups", {
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
    test_trajectories(time, state, order = order, reorder = FALSE),
    list(time, state_arr, NULL))

  ## Really simple, add an index that does not reorder anything:
  expect_equal(
    test_trajectories(time, state[1], order = order[1], reorder = TRUE),
    list(time[1], state_arr[, , , 1, drop = FALSE], NULL))
  expect_equal(
    test_trajectories(time, state[1:2], order = list(NULL, 0:6),
                      reorder = TRUE),
    list(time[1:2], state_arr[, , , 1:2, drop = FALSE], NULL))

  ## Proper reordering with the full index:
  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE),
    list(time, true, NULL))

  expect_equal(
    test_trajectories(time, state, order = order, reorder = FALSE,
                      select_particle = 3)[[2]],
    array(state_arr[, 3, , ], c(n_state, n_groups, n_time)))
  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE,
                      select_particle = 3)[[2]],
    array(true[, 3, , ], c(n_state, n_groups, n_time)))
})


test_that("can reorder trajectories on the way out", {
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
  expect_equal(test_trajectories(time, state, order = order, reorder = FALSE),
               list(time, state_arr, NULL))
  expect_equal(test_trajectories(time, state, order = order, reorder = TRUE),
               list(time, true, NULL))
})


test_that("can extract trajectories with group index, no reordering", {
  time <- seq(0, 10, length.out = 11)
  n_state <- 5
  n_particles <- 7
  n_groups <- 3
  n_time <- length(time)
  s <- lapply(seq_along(time), function(i) {
    len <- n_state * n_particles * n_groups
    # runif(n_state * n_particles * n_groups)
    array(seq_len(len) + (i - 1) * len,
          c(n_state, n_particles, n_groups))
  })

  s_arr <- array(unlist(s), c(n_state, n_particles, n_groups, n_time))

  expect_equal(test_trajectories(time, s, index_group = NULL),
               list(time, s_arr, NULL))
  expect_equal(test_trajectories(time, s, index_group = seq_len(n_groups)),
               list(time, s_arr, NULL))

  expect_equal(test_trajectories(time, s, index_group = 2),
               list(time, s_arr[, , 2, , drop = FALSE], NULL))
  expect_equal(test_trajectories(time, s, index_group = c(3, 1)),
               list(time, s_arr[, , c(3, 1), , drop = FALSE], NULL))

  m <- test_trajectories(time, s, select_particle = c(6, 4, 2))[[2]]
  expect_equal(dim(m), c(n_state, n_groups, n_time))
  expect_equal(m[, 1, ], s_arr[, 6, 1, ])
  expect_equal(m[, 2, ], s_arr[, 4, 2, ])
  expect_equal(m[, 3, ], s_arr[, 2, 3, ])

  expect_equal(
    test_trajectories(time, s,
                      index_group = c(3, 1),
                      select_particle = c(2, 6))[[2]],
    m[, c(3, 1), ])
})


test_that("can reorder trajectories on the way out", {
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
  offset <- 0
  for (i in seq_along(time)) {
    s[] <- offset + seq_along(s)
    offset <- offset + length(s)
    # s <- s + runif(length(s))
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
  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE),
    list(time, true, NULL))
  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE,
                      index_group = 1:3),
    list(time, true, NULL))

  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE,
                      index_group = 1),
    list(time, true[, , 1, , drop = FALSE], NULL))
  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE,
                      index_group = 2),
    list(time, true[, , 2, , drop = FALSE], NULL))
  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE,
                      index_group = 3),
    list(time, true[, , 3, , drop = FALSE], NULL))

  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE,
                      index_group = 3:2),
    list(time, true[, , 3:2, , drop = FALSE], NULL))

  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE,
                      index_group = 3:1),
    list(time, true[, , 3:1, ], NULL))

  m <- test_trajectories(time, state, order = order, reorder = TRUE,
                         select_particle = c(6, 4, 2))[[2]]
  expect_equal(m[, 1, ], true[, 6, 1, ])
  expect_equal(m[, 2, ], true[, 4, 2, ])
  expect_equal(m[, 3, ], true[, 2, 3, ])

  expect_equal(
    test_trajectories(time, state, order = order, reorder = TRUE,
                      index_group = c(3, 1),
                      select_particle = c(2, 6))[[2]],
    m[, c(3, 1), ])
})


## Same basic approach as the non-reordered version of trajectories.
## This is the easy bit.
test_that("can use snapshots", {
  time <- seq(0, 10, length.out = 11)
  save_snapshots <- c(3, 7)
  n_state <- 6
  n_particles <- 7
  n_groups <- 3
  n_time <- length(time)
  s <- lapply(seq_along(time), function(i) {
    array(runif(n_state * n_particles * n_groups),
          c(n_state, n_particles, n_groups))
  })
  s_arr <- array(unlist(s), c(n_state, n_particles, n_groups, n_time))
  s_arr_snapshots <- s_arr[, , , save_snapshots + 1]

  expect_equal(test_trajectories(time, s, times_snapshot = save_snapshots, reorder = FALSE),
               list(time, s_arr, s_arr_snapshots))
  expect_equal(test_trajectories(time, s, times_snapshot = save_snapshots, reorder = TRUE),
               list(time, s_arr, s_arr_snapshots))
  expect_equal(test_trajectories(time, s,
                                 order = vector("list", length(time)),
                                 times_snapshot = save_snapshots,
                                 reorder = TRUE),
               list(time, s_arr, s_arr_snapshots))

  res <- test_trajectories(time, s,
                           times_snapshot = save_snapshots,
                           select_particle = c(6, 4, 2))[[3]]
  expect_equal(dim(res), c(n_state, n_groups, 2))

  expect_equal(res[, 1, ], s_arr_snapshots[, 6, 1, ])
  expect_equal(res[, 2, ], s_arr_snapshots[, 4, 2, ])
  expect_equal(res[, 3, ], s_arr_snapshots[, 2, 3, ])
})


test_that("can reorder snapshots on the way out", {
  d <- example_trajectories(n_state = 6)

  save_snapshots <- c(3, 7)
  i <- save_snapshots + 1

  ## Pass in, but ignore index
  expect_equal(
    test_trajectories(d$time, d$state$raw,
                      times_snapshot = save_snapshots,
                      order = d$order$c,
                      reorder = FALSE),
    list(d$time, d$state$array, d$state$array[, , , i, drop = FALSE]))

  ## Still correct where we don't save trajectories
  expect_equal(
    test_trajectories(d$time, d$state$raw,
                      times_snapshot = save_snapshots,
                      order = d$order$c,
                      save_state = FALSE,
                      reorder = FALSE),
    list(numeric(), NULL, d$state$array[, , , i, drop = FALSE]))

  ## Still correct where we filter trajectories by state
  expect_equal(
    test_trajectories(d$time, d$state$raw,
                      times_snapshot = save_snapshots,
                      order = d$order$c,
                      index_state = c(2, 4, 6),
                      reorder = FALSE),
    list(d$time,
         d$state$array[c(2, 4, 6), , , , drop = FALSE],
         d$state$array[, , , i, drop = FALSE]))

  ## Reorder
  expect_equal(
    test_trajectories(d$time, d$state$raw,
                      times_snapshot = save_snapshots,
                      order = d$order$c,
                      reorder = TRUE),
    list(d$time, d$state$true, d$state$true[, , , i, drop = FALSE]))

  ## Reorder, don't save state
  expect_equal(
    test_trajectories(d$time, d$state$raw,
                      times_snapshot = save_snapshots,
                      order = d$order$c,
                      save_state = FALSE,
                      reorder = TRUE),
    list(numeric(), NULL, d$state$true[, , , i, drop = FALSE]))

  ## Reorder, filter state
  expect_equal(
    test_trajectories(d$time, d$state$raw,
                      times_snapshot = save_snapshots,
                      order = d$order$c,
                      index_state = c(2, 4, 6),
                      reorder = TRUE),
    list(d$time,
         d$state$true[c(2, 4, 6), , , , drop = FALSE],
         d$state$true[, , , i, drop = FALSE]))

  ## Select a particle
  expect_equal(
    test_trajectories(d$time, d$state$raw,
                      times_snapshot = save_snapshots,
                      order = d$order$c,
                      reorder = TRUE,
                      select_particle = 3)[[3]],
    array(d$state$true[, 3, , i], c(d$n_state, 1, 2)))
})


test_that("can reorder grouped snapshots", {
  d <- example_trajectories(n_state = 6, n_groups = 3)
  save_snapshots <- c(3, 7)
  i <- save_snapshots + 1

  expect_equal(
    test_trajectories(d$time, d$state$raw,
                      times_snapshot = save_snapshots,
                      order = d$order$c,
                      reorder = FALSE),
    list(d$time, d$state$array, d$state$array[, , , i, drop = FALSE]))

  expect_equal(
    test_trajectories(d$time, d$state$raw,
                      times_snapshot = save_snapshots,
                      order = d$order$c,
                      reorder = TRUE),
    list(d$time, d$state$true, d$state$true[, , , i, drop = FALSE]))

  res <- test_trajectories(d$time, d$state$raw,
                           times_snapshot = save_snapshots,
                           order = d$order$c,
                           select_particle = c(6, 4, 2),
                           reorder = TRUE)[[3]]
  cmp <- test_trajectories(d$time, d$state$raw,
                           times_snapshot = save_snapshots,
                           order = d$order$c,
                           reorder = TRUE)[[3]]
  expect_equal(dim(res), c(d$n_state, d$n_groups, 2))
  expect_equal(res[, 1, ], cmp[, 6, 1, ])
  expect_equal(res[, 2, ], cmp[, 4, 2, ])
  expect_equal(res[, 3, ], cmp[, 2, 3, ])
})
