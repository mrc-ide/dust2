test_that("...", {
  pars <- list(sd = 1, random_initial = TRUE)
  obj <- dust2_cpu_walk_alloc(pars, 0, 1, 10, 42, FALSE)
  ptr <- obj[[1]]
  dust2_cpu_walk_rng_state(ptr)
  dust2_cpu_walk_state(ptr)
  dust2_cpu_walk_run_steps(ptr, 3)
  dust2_cpu_walk_state(ptr)
})
