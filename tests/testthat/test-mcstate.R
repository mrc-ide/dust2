test_that("can create mcstate model from sir model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_filter_create(sir(), time_start, time, data,
                            n_particles = n_particles, seed = seed)
  packer <- mcstate2::mcstate_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  prior <- mcstate2::mcstate_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  m <- dust_filter_mcstate(obj, packer)
  expect_true(m$properties$is_stochastic)

  expect_s3_class(m, "mcstate_model")
  sampler <- mcstate2::mcstate_sampler_random_walk(diag(2) * c(0.02, 0.02))

  res <- mcstate2::mcstate_sample(m + prior, sampler, 100)
  expect_s3_class(res, "mcstate_samples")
})


test_that("can create deterministic model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  time <- c(4, 8, 12, 16)
  data <- lapply(1:4, function(i) list(incidence = i))
  dt <- 1
  n_particles <- 100
  seed <- 42

  obj <- dust_unfilter_create(sir(), time_start, time, data)
  packer <- mcstate2::mcstate_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  prior <- mcstate2::mcstate_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  m <- dust_filter_mcstate(obj, packer)
  expect_false(m$properties$is_stochastic)

  expect_s3_class(m, "mcstate_model")
  sampler <- mcstate2::mcstate_sampler_random_walk(diag(2) * c(0.02, 0.02))
  res <- mcstate2::mcstate_sample(m + prior, sampler, 10)
  expect_s3_class(res, "mcstate_samples")
})