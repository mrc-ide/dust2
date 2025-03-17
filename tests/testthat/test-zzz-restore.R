test_that("can serialise and restore model", {
  mysir <- local_sir_generator()

  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  set.seed(1)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  obj <- dust_filter_create(mysir, time_start, data, n_particles = 100,
                            seed = 42)
  m <- dust_likelihood_monty(obj, packer) + prior
  s <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))

  set.seed(1)
  res1 <- monty::monty_sample(m, s, 20, initial = c(0.2, 0.1))

  set.seed(1)
  path <- withr::local_tempdir()
  monty::monty_sample_manual_prepare(m, s, 20, path, initial = c(0.2, 0.1))
  monty::monty_sample_manual_run(1, path)
  res2 <- monty::monty_sample_manual_collect(path)
  expect_equal(res2, res1)
})


test_that("Can restore model in monty context", {
  mysir <- local_sir_generator()

  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  set.seed(1)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)

  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  obj <- dust_filter_create(mysir, time_start, data, n_particles = 100)
  m <- dust_likelihood_monty(obj, packer) + prior
  s <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))
  r <- monty::monty_runner_callr(2)

  set.seed(1)
  res <- monty::monty_sample(m, s, 20, initial = c(0.2, 0.1), runner = r)

  set.seed(1)
  cmp <- monty::monty_sample(m, s, 20, initial = c(0.2, 0.1))

  expect_equal(res, cmp)
})
