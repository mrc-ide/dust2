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


test_that("can get trajectories from model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                            seed = 42)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  m <- dust_likelihood_monty(obj, packer, save_trajectories = TRUE)
  expect_true(m$properties$has_observer)

  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  posterior <- m + prior

  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))
  res <- monty::monty_sample(posterior, sampler, 27, initial = c(.2, .1),
                             n_chains = 3)

  expect_equal(names(res$observations), "trajectories")
  expect_equal(dim(res$observations$trajectories), c(5, 4, 27, 3))
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
  m1 <- dust_likelihood_monty(obj, packer,
                              save_trajectories = c("I", "cases_inc"))
  p1 <- m1 + prior
  res1 <- monty::monty_sample(p1, sampler, 13, initial = c(.2, .1),
                              n_chains = 3)

  set.seed(1)
  m2 <- dust_likelihood_monty(obj, packer, save_trajectories = TRUE)
  p2 <- m2 + prior
  res2 <- monty::monty_sample(p2, sampler, 13, initial = c(.2, .1),
                              n_chains = 3)

  expect_equal(res1$pars, res2$pars)
  expect_equal(dim(res1$observations$trajectories), c(2, 4, 13, 3))
  expect_equal(res1$observations$trajectories,
               res2$observations$trajectories[c(2, 5), , , ])
  expect_equal(names(res1$observations), "trajectories")
})


test_that("can subset trajectories from deterministic model", {
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

  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))

  set.seed(1)
  m1 <- dust_likelihood_monty(obj, packer,
                              save_trajectories = c("I", "cases_inc"))
  p1 <- m1 + prior
  res1 <- monty::monty_sample(p1, sampler, 13, initial = c(.2, .1),
                              n_chains = 3)

  set.seed(1)
  m2 <- dust_likelihood_monty(obj, packer, save_trajectories = TRUE)
  p2 <- m2 + prior
  res2 <- monty::monty_sample(p2, sampler, 13, initial = c(.2, .1),
                              n_chains = 3)

  expect_equal(dim(res1$observations$trajectories), c(2, 4, 13, 3))
  expect_equal(res1$observations$trajectories,
               res2$observations$trajectories[c(2, 5), , , ])
  expect_equal(names(res1$observations), "trajectories")
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
                             save_trajectories = TRUE)
  expect_true(m$properties$has_observer)
  res <- monty::monty_sample(m + prior, sampler, 13, n_chains = 3)
  expect_equal(names(res$observations), c("state", "trajectories"))
  expect_equal(dim(res$observations$state), c(5, 13, 3))
  expect_equal(dim(res$observations$trajectories), c(5, 4, 13, 3))
  expect_equal(res$observations$state, res$observations$trajectories[, 4, , ])
})


test_that("validate request to save trajectories", {
  expect_equal(validate_save_trajectories(FALSE), list(enabled = FALSE))
  expect_equal(validate_save_trajectories(TRUE),
               list(enabled = TRUE, index = NULL, subset = NULL))
  expect_error(validate_save_trajectories(NULL),
               "Invalid value for 'save_trajectories'")
  expect_error(validate_save_trajectories(1:3),
               "Invalid value for 'save_trajectories'")
  expect_equal(validate_save_trajectories(c("a", "b")),
               list(enabled = TRUE, index = NULL, subset = c("a", "b")))
})


test_that("can use names for groups", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)

  time_start <- 0
  data1 <- data.frame(group = rep(1:2, each = 4),
                      time = rep(c(4, 8, 12, 16), 2),
                      incidence = c(1:4, 2:5))
  data2 <- data1
  data2$group <- letters[data2$group]

  obj1 <- dust_filter_create(sir(), time_start, data1, n_particles = 100,
                             seed = 42)
  obj2 <- dust_filter_create(sir(), time_start, data2, n_particles = 100,
                             seed = 42, n_groups = 2)

  packer <- monty::monty_packer_grouped(
    c("a", "b"),
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))

  m2 <- dust_likelihood_monty(obj2, packer)
  p2 <- c(0.2, 0.1, 0.25, 0.1)

  ll1 <- dust_likelihood_run(obj1, unname(packer$unpack(p2)))
  ll2 <- monty::monty_model_density(m2, p2)
  expect_equal(ll2, sum(ll1))
})


test_that("can run simple grouped unfilter", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  
  time_start <- 0
  data1 <- data.frame(group = rep(1:2, each = 4),
                      time = rep(c(4, 8, 12, 16), 2),
                      incidence = c(1:4, 2:5))
  data2 <- data1
  data2$group <- letters[data2$group]
  
  obj1 <- dust_unfilter_create(sir(), time_start, data1)
  obj2 <- dust_unfilter_create(sir(), time_start, data2, n_groups = 2)
  
  packer <- monty::monty_packer_grouped(
    c("a", "b"),
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  
  m2 <- dust_likelihood_monty(obj2, packer)
  p2 <- c(0.2, 0.1, 0.25, 0.1)
  
  ll1 <- dust_likelihood_run(obj1, unname(packer$unpack(p2)))
  ll2 <- monty::monty_model_density(m2, p2)
  expect_equal(ll2, sum(ll1))
})


test_that("observers do not change the steps", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj1 <- dust_filter_create(sir(), time_start, data, n_particles = 100)
  obj2 <- dust_filter_create(sir(), time_start, data, n_particles = 100)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  m1 <- dust_likelihood_monty(obj1, packer, save_trajectories = TRUE)
  m2 <- dust_likelihood_monty(obj2, packer, save_trajectories = FALSE)

  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  posterior1 <- m1 + prior
  posterior2 <- m2 + prior

  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))
  set.seed(1)
  res1 <- monty::monty_sample(posterior1, sampler, 27, initial = c(.2, .1),
                              n_chains = 3)
  set.seed(1)
  res2 <- monty::monty_sample(posterior2, sampler, 27, initial = c(.2, .1),
                              n_chains = 3)
  expect_equal(res1$pars, res2$pars)
  expect_equal(dim(res1$observations$trajectories), c(5, 4, 27, 3))
  expect_null(res2$observations$trajectories)
})


test_that("can get snapshots from model", {
  pars <- list(beta = 0.1, gamma = 0.2, N = 1000, I0 = 10, exp_noise = 1e6)
  time_start <- 0
  data <- data.frame(time = c(4, 8, 12, 16), incidence = 1:4)
  obj <- dust_filter_create(sir(), time_start, data, n_particles = 100,
                            seed = 42)
  packer <- monty::monty_packer(
    c("beta", "gamma"),
    fixed = list(N = 1000, I0 = 10, exp_noise = 1e6))
  m <- dust_likelihood_monty(obj, packer, save_snapshots = c(4, 12))
  expect_true(m$properties$has_observer)

  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 0.5)
    gamma ~ Exponential(mean = 0.5)
  })

  posterior <- m + prior

  sampler <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.02))
  res <- monty::monty_sample(posterior, sampler, 27, initial = c(.2, .1),
                             n_chains = 3)

  expect_equal(names(res$observations), "snapshots")
  expect_equal(dim(res$observations$snapshots), c(5, 2, 27, 3))
})


test_that("cope with changing the trajectory index", {
  d <- data.frame(
    time = 1:5,
    incidence = c(12, 23, 25, 36, 30))

  filter <- dust2::dust_filter_create(sir, time_start = 0,
                                      data = d, dt = 0.25,
                                      n_particles = 200)

  packer <- monty::monty_packer(scalar = c("beta", "gamma"),
                                fixed = list(I0 = 10, N = 1000))
  vcv <- diag(2) * 0.01
  sampler <- monty::monty_sampler_random_walk(vcv)
  prior <- monty::monty_dsl({
    beta ~ Exponential(mean = 1)
    gamma ~ Exponential(mean = 0.5)
  })

  likelihood <- dust2::dust_likelihood_monty(filter, packer,
                                             save_trajectories = TRUE)
  posterior <- likelihood + prior
  samples <- monty::monty_sample(posterior, sampler, n_steps = 3,
                                 initial = c(0.3, 0.1),
                                 n_chains = 2)
  expect_equal(
    dim(samples$observations$trajectories),
    c(5, 5, 3, 2)) # state, time, steps, chains

  likelihood <- dust2::dust_likelihood_monty(filter, packer,
                                             save_trajectories = "cases_inc")
  posterior <- likelihood + prior
  samples <- monty::monty_sample(posterior, sampler, n_steps = 3,
                                 initial = c(0.3, 0.1),
                                 n_chains = 2)
  expect_equal(
    dim(samples$observations$trajectories),
    c(1, 5, 3, 2)) # state, time, steps, chains
})
