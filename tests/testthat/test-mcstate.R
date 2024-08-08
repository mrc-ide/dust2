test_that("can create mcstate model from sir model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                            seed = 42)
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
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sir(), time_start, data)
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


test_that("can avoid errors by converting to impossible density", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  obj <- dust_unfilter_create(sir(), time_start, data)
  packer <- mcstate2::mcstate_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  prior <- mcstate2::mcstate_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  m1 <- dust_filter_mcstate(obj, packer)
  m2 <- dust_filter_mcstate(obj, packer, failure_is_impossible = TRUE)

  expect_error(m1$density(c(-1, -1)),
               "Invalid call to binomial")
  expect_equal(m2$density(c(-1, -1)), -Inf)
})
