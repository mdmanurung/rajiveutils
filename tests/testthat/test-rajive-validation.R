# Input validation and warning tests for Rajive (W-M10, W-M11).

test_that("Rajive aborts on near-constant columns", {
  set.seed(12)
  n <- 30
  X1 <- matrix(rnorm(n * 8), n, 8)
  X2 <- matrix(rnorm(n * 7), n, 7)

  # Force degeneracy in one feature.
  X1[, 1] <- 5

  expect_error(
    Rajive(list(X1, X2), initial_signal_ranks = c(3, 3)),
    class = "rajiveplus_degenerate_block"
  )
})

test_that("Rajive warns when n < sum(initial_signal_ranks)", {
  set.seed(13)
  n <- 12
  X1 <- matrix(rnorm(n * 10), n, 10)
  X2 <- matrix(rnorm(n * 9), n, 9)

  expect_warning(
    Rajive(list(X1, X2), initial_signal_ranks = c(7, 6), joint_rank = 1,
           n_wedin_samples = 20, n_rand_dir_samples = 20, num_cores = 1),
    class = "rajiveplus_underdetermined"
  )
})

test_that("Rajive adjusts RNG kind for parallel and restores after exit", {
  skip_on_cran()

  set.seed(14)
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)

  # Force a non-L'Ecuyer kind where available.
  RNGkind("Mersenne-Twister")

  X1 <- matrix(rnorm(25 * 8), 25, 8)
  X2 <- matrix(rnorm(25 * 7), 25, 7)

  expect_warning(
    Rajive(list(X1, X2), initial_signal_ranks = c(3, 3), joint_rank = 1,
           n_wedin_samples = 20, n_rand_dir_samples = 20, num_cores = 2),
    class = "rajiveplus_rngkind_adjusted"
  )

  # RNG kind should be restored on exit.
  expect_identical(RNGkind(), c("Mersenne-Twister", old_kind[2], old_kind[3]))
})
