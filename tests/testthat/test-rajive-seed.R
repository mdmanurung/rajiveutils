# Tests for W-R7: seed argument on Rajive().

test_that("Rajive with same seed produces bit-identical results", {
  skip_on_cran()

  n   <- 30
  set.seed(0)
  X1  <- matrix(rnorm(n * 8), n, 8)
  X2  <- matrix(rnorm(n * 7), n, 7)

  fit1 <- Rajive(
    list(X1, X2),
    initial_signal_ranks = c(3, 3),
    n_wedin_samples      = 50,
    n_rand_dir_samples   = 50,
    num_cores            = 1L,
    seed                 = 42L
  )

  fit2 <- Rajive(
    list(X1, X2),
    initial_signal_ranks = c(3, 3),
    n_wedin_samples      = 50,
    n_rand_dir_samples   = 50,
    num_cores            = 1L,
    seed                 = 42L
  )

  # Joint rank and scores must be identical.
  expect_equal(fit1$joint_rank, fit2$joint_rank)
  expect_equal(fit1$joint_scores, fit2$joint_scores)

  # Wedin samples must be identical.
  ws1 <- fit1$joint_rank_sel[["wedin"]]$wedin_samples
  ws2 <- fit2$joint_rank_sel[["wedin"]]$wedin_samples
  expect_equal(ws1, ws2)
})

test_that("Rajive with different seeds returns valid bound samples", {
  skip_on_cran()

  n   <- 30
  set.seed(1)
  X1  <- matrix(rnorm(n * 8), n, 8)
  X2  <- matrix(rnorm(n * 7), n, 7)

  fit1 <- Rajive(
    list(X1, X2),
    initial_signal_ranks = c(3, 3),
    n_wedin_samples      = 50,
    n_rand_dir_samples   = 50,
    num_cores            = 1L,
    seed                 = 1L
  )

  fit2 <- Rajive(
    list(X1, X2),
    initial_signal_ranks = c(3, 3),
    n_wedin_samples      = 50,
    n_rand_dir_samples   = 50,
    num_cores            = 1L,
    seed                 = 2L
  )

  ws1 <- fit1$joint_rank_sel[["wedin"]]$wedin_samples
  ws2 <- fit2$joint_rank_sel[["wedin"]]$wedin_samples
  expect_equal(length(ws1), length(ws2))
  expect_true(all(is.finite(ws1)))
  expect_true(all(is.finite(ws2)))
})
