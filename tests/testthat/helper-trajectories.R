example_trajectories <- function(time = 0:10, n_state = 3, n_particles = 7,
                                 n_groups = 1) {
  time <- seq(0, 10, length.out = 11)
  n_time <- length(time)

  state <- vector("list", length(time))
  order <- vector("list", length(time))
  true <- array(NA_real_, c(n_state, n_particles, n_groups, n_time))
  s <- array(0, c(n_state, n_particles, n_groups))
  set.seed(1)
  offset <- 0
  for (i in seq_along(time)) {
    s[] <- offset + seq_along(s)
    offset <- offset + length(s)
    s <- s + runif(length(s))
    if (i > 1 && i %% 2 == 0) {
      k <- replicate(n_groups, sample(n_particles, replace = TRUE))
      order[[i]] <- as.integer(k)
      for (j in seq_len(n_groups)) {
        s[, , j] <- s[, k[, j], j] + 1
        true[, , j, seq_len(i - 1)] <- true[, k[, j], j, seq_len(i - 1)]
      }
    }
    state[[i]] <- s
    true[, , , i] <- s
  }

  state_arr <- array(unlist(state), c(n_state, n_particles, n_groups, n_time))

  order_c <- lapply(order, function(x) if (!is.null(x)) x - 1L)

  list(
    time = time,
    state = list(true = true, array = state_arr, raw = state),
    reorder = !vlapply(order, is.null),
    order = list(c = order_c, r = order),
    n_time = n_time,
    n_state = n_state,
    n_particles = n_particles,
    n_groups = n_groups)
}


## For testing, these are reordering functions that will work with the
## output of the above function to reorder state or snapshots for a
## single or all particles.
cmp_reorder_state_single <- function(idx_particle, d) {
  ret <- array(NA_real_, c(d$n_state, d$n_time))
  for (i in rev(seq_len(d$n_time))) {
    ret[, i] <- d$state$raw[[i]][, idx_particle, 1]
    if (d$reorder[[i]]) {
      idx_particle <- d$order$r[[i]][idx_particle]
    }
  }
  ret
}


cmp_reorder_state_all <- function(d) {
  idx_particle <- seq_len(d$n_particles)
  ret <- array(NA_real_, c(d$n_state, d$n_particles, d$n_time))
  for (i in rev(seq_len(d$n_time))) {
    ret[, , i] <- d$state$raw[[i]][, idx_particle, 1]
    if (d$reorder[[i]]) {
      idx_particle <- d$order$r[[i]][idx_particle]
    }
  }
  ret
}


cmp_reorder_snapshot_single <- function(idx_particle, t_snapshot, d) {
  d_snapshot <- d$state$raw[match(t_snapshot, d$time)]
  n_snapshots <- length(t_snapshot)

  ret <- array(NA_real_, c(d$n_state, n_snapshots))
  it_t <- n_snapshots
  for (i in rev(seq_len(d$n_time))) {
    if (length(idx_particle) != 1) browser()
    t <- d$time[[i]]
    if (t == t_snapshot[[it_t]]) {
      ret[, it_t] <- d_snapshot[[it_t]][, idx_particle, 1L]
      if (it_t == 1) {
        break
      }
      it_t <- it_t - 1L
    }

    if (d$reorder[[i]]) {
      idx_particle <- d$order$r[[i]][idx_particle]
    }
  }
  ret
}


cmp_reorder_snapshot_all <- function(t_snapshot, d) {
  d_snapshot <- d$state$raw[match(t_snapshot, d$time)]
  n_snapshots <- length(t_snapshot)

  idx_particle <- seq_len(d$n_particles)

  ret <- array(NA_real_, c(d$n_state, d$n_particles, n_snapshots))
  ## cmp <- array(NA_real_, c(d$n_state, d$n_particles, d$n_time))
  it_t <- n_snapshots

  for (i in rev(seq_len(d$n_time))) {
    t <- d$time[[i]]
    if (t == t_snapshot[[it_t]]) {
      ret[, , it_t] <- d_snapshot[[it_t]][, idx_particle, 1L]
      if (it_t == 1) {
        break
      }
      it_t <- it_t - 1L
    }

    if (d$reorder[[i]]) {
      idx_particle <- d$order$r[[i]][idx_particle]
    }
  }
  ret
}
