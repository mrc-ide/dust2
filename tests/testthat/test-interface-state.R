## There are many tests for checking that state conforms to a system,
## so they're all here so that we can have a nice set of labelled
## tests.  This is quite repetitive (as is the implementation) but the
## aim is to convey some hint to the user about what they have done
## wrong.

test_that("can provide a vector to set into a vector system", {
  expect_equal(
    prepare_state(1:3, NULL, NULL, NULL, 3, 1, 1, FALSE, FALSE),
    list(state = 1:3,
         index_state = NULL,
         index_particle = NULL,
         index_group = NULL,
         recycle_particle = FALSE,
         recycle_group = FALSE))
})

test_that("can provide a vector with index to set into vector system", {
  expect_equal(
    prepare_state(1:3, 3:5, NULL, NULL, 5, 1, 1, FALSE, FALSE),
    list(state = 1:3,
         index_state = 3:5,
         index_particle = NULL,
         index_group = NULL,
         recycle_particle = FALSE,
         recycle_group = FALSE))
})


test_that("error if incorrect data provided to vector system", {
  expect_error(
    prepare_state(cbind(1:3), NULL, NULL, NULL, 3, 1, 1, FALSE, FALSE),
    "Expected 'state' to be a vector but was given a matrix")
  expect_error(
    prepare_state(array(1:3, c(3, 1, 1)), NULL, NULL, NULL, 3, 1, 1,
                  FALSE, FALSE),
    "Expected 'state' to be a vector but was given a 3-dimensional array")
})


test_that("validate that the index provided for state is reasonable", {
  expect_error(
    prepare_state(1:4, 3:6, NULL, NULL, 5, 1, 1, FALSE, FALSE, name = "state"),
    "All elements of 'index_state' must be at most 5")
})


## TODO: hint here about if we have provided an index, and what we
## did recieve.
test_that("error if incorrect length data provided to vector system", {
  expect_error(
    prepare_state(1:3, NULL, NULL, NULL, 5, 1, 1, FALSE, FALSE, name = "state"),
    "Expected 'state' to have length 5")
  expect_error(
    prepare_state(1:3, 3:6, NULL, NULL, 6, 1, 1, FALSE, FALSE, name = "state"),
    "Expected 'state' to have length 4")
})


test_that("can provide a vector to set into matrix (s x p) system", {
  expect_equal(
    prepare_state(1:3, NULL, NULL, NULL, 3, 5, 1, TRUE, FALSE),
    list(state = 1:3,
         index_state = NULL,
         index_particle = NULL,
         index_group = NULL,
         recycle_particle = TRUE,
         recycle_group = FALSE))
})


test_that("can provide a matrix to set into a matrix (s x p) sytem", {
  expect_equal(
    prepare_state(cbind(1:3), NULL, NULL, NULL, 3, 5, 1, TRUE, FALSE),
    list(state = cbind(1:3),
         index_state = NULL,
         index_particle = NULL,
         index_group = NULL,
         recycle_particle = TRUE,
         recycle_group = FALSE))
  state <- matrix(1:15, 3)
  expect_equal(
    prepare_state(matrix(1:15, 3), NULL, NULL, NULL, 3, 5, 1, TRUE, FALSE),
        list(state = state,
         index_state = NULL,
         index_particle = NULL,
         index_group = NULL,
         recycle_particle = FALSE,
         recycle_group = FALSE))
})


test_that("can provide state with state index into matrix system", {
  expect_equal(
    prepare_state(1:3, 3:5, NULL, NULL, 7, 5, 1, TRUE, FALSE),
    list(state = 1:3,
         index_state = 3:5,
         index_particle = NULL,
         index_group = NULL,
         recycle_particle = TRUE,
         recycle_group = FALSE))
  state <- matrix(1:15, 3:5)
  expect_equal(
    prepare_state(state, 3:5, NULL, NULL, 7, 5, 1, TRUE, FALSE),
    list(state = state,
         index_state = 3:5,
         index_particle = NULL,
         index_group = NULL,
         recycle_particle = FALSE,
         recycle_group = FALSE))
})


test_that("can provide state with particle index into matrix system", {
  expect_equal(
    prepare_state(1:3, NULL, 2:3, NULL, 3, 5, 1, TRUE, FALSE),
    list(state = 1:3,
         index_state = NULL,
         index_particle = 2:3,
         index_group = NULL,
         recycle_particle = TRUE,
         recycle_group = FALSE))
  state <- matrix(1:6, 3)
  expect_equal(
    prepare_state(state, NULL, 2:3, NULL, 3, 5, 1, TRUE, FALSE),
    list(state = state,
         index_state = NULL,
         index_particle = 2:3,
         index_group = NULL,
         recycle_particle = FALSE,
         recycle_group = FALSE))
})


test_that("validate that the index provided for particle is reasonable", {
  expect_error(
    prepare_state(1:4, NULL, 5:6, NULL, 4, 2, 1, TRUE, FALSE),
    "All elements of 'index_particle' must be at most 2")
  expect_error(
    prepare_state(1:4, NULL, c(-1, 0), NULL, 4, 2, 1, TRUE, FALSE),
    "All elements of 'index_particle' must be at least 1")
  expect_error(
    prepare_state(1:4, NULL, c(1, 1, 2), NULL, 4, 2, 1, TRUE, FALSE),
    "All elements of 'index_particle' must be distinct")
})


test_that("can set state from vector", {
  sys <- dust_system_create(sir(), list(), n_particles = 3)
  dust_system_set_state(sys, as.numeric(1:5))
  expect_equal(dust_system_state(sys), matrix(1:5, 5, 3))
})


test_that("can set state from matrix", {
  sys <- dust_system_create(sir(), list(), n_particles = 3)
  m <- matrix(runif(15), 5, 3)
  dust_system_set_state(sys, m)
  expect_equal(dust_system_state(sys), m)
})


test_that("can set a fraction of states", {
  sys <- dust_system_create(sir(), list(), n_particles = 3)
  m <- matrix(runif(15), 5, 3)
  dust_system_set_state(sys, m)
  m2 <- matrix(seq(16, length.out = 6), 2, 3)
  dust_system_set_state(sys, m2, index_state = c(2, 4))
  m[c(2, 4), ] <- m2
  expect_equal(dust_system_state(sys), m)
})


test_that("can set a fraction of states from a vector", {
  sys <- dust_system_create(sir(), list(), n_particles = 3)
  m <- matrix(runif(15), 5, 3)
  dust_system_set_state(sys, m)
  m2 <- c(16, 17)
  dust_system_set_state(sys, m2, index_state = c(2, 4))
  m[c(2, 4), ] <- m2
  expect_equal(dust_system_state(sys), m)
})


test_that("can set a fraction of states from a scalar", {
  sys <- dust_system_create(sir(), list(), n_particles = 3)
  m <- matrix(runif(15), 5, 3)
  dust_system_set_state(sys, m)
  m2 <- 16
  dust_system_set_state(sys, m2, index_state = 2)
  m[2, ] <- m2
  expect_equal(dust_system_state(sys), m)
})


test_that("can set a fraction of particles", {
  sys <- dust_system_create(sir(), list(), n_particles = 6)
  m <- matrix(as.numeric(1:30), 5, 6)
  dust_system_set_state(sys, m)
  m2 <- matrix(seq(31, length.out = 15), 5, 3)
  dust_system_set_state(sys, m2, index_particle = c(2, 4, 6))
  m[, c(2, 4, 6)] <- m2
  expect_equal(dust_system_state(sys), m)
})


test_that("can set a fraction of particles from a vector", {
  sys <- dust_system_create(sir(), list(), n_particles = 6)
  m <- matrix(as.numeric(1:30), 5, 6)
  dust_system_set_state(sys, m)
  m2 <- seq(31, length.out = 5)
  dust_system_set_state(sys, m2, index_particle = c(2, 4, 6))
  m[, c(2, 4, 6)] <- m2
  expect_equal(dust_system_state(sys), m)
})


test_that("can set a fraction of both particles and states", {
  sys <- dust_system_create(sir(), list(), n_particles = 6)
  m <- matrix(as.numeric(1:30), 5, 6)
  dust_system_set_state(sys, m)
  m2 <- matrix(seq(31, length.out = 6), 2, 3)
  dust_system_set_state(sys, m2,
                        index_state = c(2, 4),
                        index_particle = c(2, 4, 6))
  m[c(2, 4), c(2, 4, 6)] <- m2
  expect_equal(dust_system_state(sys), m)
})


test_that("can set vector into grouped sytem", {
  sys <- dust_system_create(walk(), rep(list(list(len = 5, sd = 1)), 3),
                            n_particles = 4,
                            n_groups = 3)
  dust_system_set_state(sys, 1:5)
  expect_equal(dust_system_state(sys), array(1:5, c(5, 4, 3)))
})


test_that("can set matrix into grouped sytem", {
  sys <- dust_system_create(walk(), rep(list(list(len = 5, sd = 1)), 3),
                            n_particles = 4,
                            n_groups = 3)
  m <- matrix(1:20, 5, 4)
  dust_system_set_state(sys, m)
  expect_equal(dust_system_state(sys), array(m, c(5, 4, 3)))
})


test_that("can set array into grouped sytem", {
  sys <- dust_system_create(walk(), rep(list(list(len = 5, sd = 1)), 3),
                            n_particles = 4,
                            n_groups = 3)
  m <- array(1:60, c(5, 4, 3))
  dust_system_set_state(sys, m)
  expect_equal(dust_system_state(sys), m)
})


test_that("can set subset into grouped system", {
  sys <- dust_system_create(walk(), rep(list(list(len = 5, sd = 1)), 3),
                            n_particles = 4,
                            n_groups = 3)
  m1 <- array(1:60, c(5, 4, 3))
  dust_system_set_state(sys, m1)

  i <- c(2, 4)
  j <- c(1, 3)
  k <- c(2, 3)
  m2 <- array(100 + seq_len(8), c(2, 2, 2))
  dust_system_set_state(sys, m2, index_state = i, index_particle = j,
                        index_group = k)
  m1[i, j, k] <- m2
  expect_equal(dust_system_state(sys), m1)
})
