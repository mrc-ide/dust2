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
})


test_that("can use history", {
  time <- seq(0, 10, length.out = 11)
  n_state <- 6
  n_particles <- 7
  n_groups <- 3
  s <- lapply(seq_along(time), function(i) {
    array(runif(n_state * n_particles * n_groups),
          c(n_state, n_particles, n_groups))
  })
  expect_equal(test_history(time, s, NULL, TRUE),
               list(time, unlist(s)))
  expect_equal(test_history(time, s, NULL, FALSE),
               list(time, unlist(s)))
  expect_equal(test_history(time, s[1:3], NULL, TRUE),
               list(time[1:3], unlist(s[1:3])))
  expect_equal(test_history(time, s, vector("list", length(time)), TRUE),
               list(time, unlist(s)))



})


test_that("can reorder history on the way out", {
  ## This is really hard to get right so let's actually simulate forward:
  time <- seq(0, 10, length.out = 11)
  n_state <- 6
  n_particles <- 7
  n_groups <- 3
  state <- vector("list", length(time))
  order <- vector("list", length(time))
  for (i in seq_along(time)) {
    if (i == 1) {
      s <- array(runif(n_state * n_particles * n_groups),
                 c(n_state, n_particles, n_groups))
    } else if (i %% 2 == 0) {
      s <- s + 1
    } else {
      k <- replicate(n_groups, sample(n_particles, replace = TRUE))
      order[[i]] <- as.integer(k - 1L)
      for (j in seq_len(n_groups)) {
        s[, , j] <- s[, k[, j], j] + 1
      }
    }
    state[[i]] <- s
  }

  expect_equal(test_history(time, state, order, FALSE),
               list(time, unlist(state)))
})
