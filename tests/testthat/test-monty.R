test_that("can create monty model from sir model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                            seed = 42)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  m <- dust_likelihood_monty(obj, packer)
  expect_true(m$properties$is_stochastic)

  expect_s3_class(m, "monty_model")
  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))

  res <- monty::monty_sample(m + prior, sampler, 100)
  expect_s3_class(res, "monty_samples")
})


test_that("can create deterministic model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sir(), time_start, data)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  m <- dust_likelihood_monty(obj, packer)
  expect_false(m$properties$is_stochastic)
  expect_true(m$properties$has_gradient)

  expect_length(m$density(c(0.2, 0.1)), 1)
  expect_length(m$gradient(c(0.2, 0.1)), 3)

  expect_s3_class(m, "monty_model")
  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))
  res <- monty::monty_sample(m + prior, sampler, 10)
  expect_s3_class(res, "monty_samples")
})


test_that("can avoid errors by converting to impossible density", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sir(), time_start, data)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  m1 <- dust_likelihood_monty(obj, packer)
  m2 <- dust_likelihood_monty(obj, packer, failure_is_impossible = TRUE)

  expect_error(m1$density(c(-1, -1)),
               "Invalid call to binomial")
  expect_equal(m2$density(c(-1, -1)), -Inf)
})


test_that("can create wrapper around filter with multiple pars", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(group = rep(1:2, each = 4),
                     time = rep(c(4, 8, 12, 16), 2),
                     incidence = c(1:4, 2:5))

  obj <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                            seed = 42)
  obj1 <- dust_filter_create(sir(), time_start, data[data$group == 1, ],
                             n_particles = 100, seed = 42)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))

  m <- dust_likelihood_monty(obj, packer)
  expect_true(m$properties$allow_multiple_parameters)
  expect_true(m$properties$is_stochastic)

  m1 <- dust_likelihood_monty(obj1, packer)
  expect_false(m1$properties$allow_multiple_parameters)

  p <- cbind(c(0.2, 0.1), c(0.25, 0.1))

  ll <- monty::monty_model_density(m, p)
  ll1 <- monty::monty_model_density(m1, p[, 1])
  expect_length(ll, 2)
  expect_equal(ll[[1]], ll1)
})


test_that("can get trajectories from model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                            seed = 42)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  m <- dust_likelihood_monty(obj, packer, save_history = TRUE)
  expect_true(m$properties$has_observer)

  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  posterior <- m + prior

  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))
  res <- monty::monty_sample(posterior, sampler, 27, initial = c(.2, .1),
                             n_chains = 3)

  expect_equal(names(res$observations), "history")
  expect_equal(dim(res$observations$history), c(5, 4, 27, 3))
})


test_that("can subset trajectories from model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                            seed = 42)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))

  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))

  set.seed(1)
  m1 <- dust_likelihood_monty(obj, packer, save_history = c("I", "cases_inc"))
  p1 <- m1 + prior
  res1 <- monty::monty_sample(p1, sampler, 13, initial = c(.2, .1),
                              n_chains = 3)

  set.seed(1)
  m2 <- dust_likelihood_monty(obj, packer, save_history = TRUE)
  p2 <- m2 + prior
  res2 <- monty::monty_sample(p2, sampler, 13, initial = c(.2, .1),
                              n_chains = 3)

  expect_equal(dim(res1$observations$history), c(2, 4, 13, 3))
  expect_equal(res1$observations$history,
               res2$observations$history[c(2, 5), , , ])
  expect_equal(names(res1$observations), "history")
})


test_that("can record final state", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                            seed = 42)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))

  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))

  set.seed(1)
  m <- dust_likelihood_monty(obj, packer, save_state = TRUE)
  expect_true(m$properties$has_observer)
  res <- monty::monty_sample(m + prior, sampler, 13, n_chains = 3)
  expect_equal(names(res$observations), "state")
  expect_equal(dim(res$observations$state), c(5, 13, 3))
})


test_that("can record both state and trajectories", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                            seed = 42)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))

  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))

  set.seed(1)
  m <- dust_likelihood_monty(obj, packer, save_state = TRUE,
                             save_history = TRUE)
  expect_true(m$properties$has_observer)
  res <- monty::monty_sample(m + prior, sampler, 13, n_chains = 3)
  expect_equal(names(res$observations), c("state", "history"))
  expect_equal(dim(res$observations$state), c(5, 13, 3))
  expect_equal(dim(res$observations$history), c(5, 4, 13, 3))
  expect_equal(res$observations$state, res$observations$history[, 4, , ])
})
