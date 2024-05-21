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
