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

  m <- dust_filter_monty(obj, packer)
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

  m <- dust_filter_monty(obj, packer)
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

  m1 <- dust_filter_monty(obj, packer)
  m2 <- dust_filter_monty(obj, packer, failure_is_impossible = TRUE)

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

  m <- dust_filter_monty(obj, packer)
  expect_true(m$properties$allow_multiple_parameters)
  expect_true(m$properties$is_stochastic)

  m1 <- dust_filter_monty(obj1, packer)
  expect_false(m1$properties$allow_multiple_parameters)

  p <- cbind(c(0.2, 0.1), c(0.25, 0.1))

  ll <- monty::monty_model_density(m, p)
  ll1 <- monty::monty_model_density(m1, p[, 1])
  expect_length(ll, 2)
  expect_equal(ll[[1]], ll1)
})
