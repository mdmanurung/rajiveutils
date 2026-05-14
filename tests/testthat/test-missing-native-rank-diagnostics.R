test_that("candidate-rank diagnostics report observed prediction and support metrics", {
  set.seed(8701)
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[2, 3] <- FALSE

  diag <- diagnose_missing_ranks(blocks, candidates = 0:2,
                                 initial_signal_ranks = c(2L, 2L),
                                 mask = mask)

  expect_equal(names(diag), c("joint_rank", "prediction_error",
                              "weak_support_rate", "missing_fraction",
                              "composite_score"))
  expect_equal(diag$missing_fraction, rep(mean(!unlist(mask)), 3))
  expect_true(all(diag$prediction_error >= 0))
})

test_that("rank diagnostics report weak support rate", {
  set.seed(8702)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block2[4, ] <- FALSE

  diag <- diagnose_missing_ranks(blocks, candidates = c(1L),
                                 initial_signal_ranks = c(2L, 2L),
                                 mask = mask)
  expect_gt(diag$weak_support_rate, 0)
})

test_that("rank diagnostics report prediction error on the fitted preprocessing scale", {
  set.seed(8703)
  blocks <- list(
    block1 = matrix(rnorm(72, mean = 100, sd = 20), nrow = 12),
    block2 = matrix(rnorm(60, mean = -75, sd = 10), nrow = 12)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[c(2, 7), 3] <- FALSE
  mask$block2[4, 2] <- FALSE
  control <- rajive_missing_control(center = TRUE, scale = TRUE,
                                    normalize = TRUE)

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L, seed = 8703,
                missing_control = control,
                identifiability_norm = "l2")
  diag <- diagnose_missing_ranks(blocks, candidates = 1L,
                                 initial_signal_ranks = c(2L, 2L),
                                 mask = mask,
                                 missing_control = control,
                                 seed = 8703,
                                 identifiability_norm = "l2")

  prep <- rajiveplus:::.center_scale_observed(
    blocks, mask, center = TRUE, scale = TRUE, normalize = TRUE
  )
  observed_ss <- sum(vapply(seq_along(prep$blocks), function(k) {
    sum(prep$blocks[[k]][mask[[k]]]^2)
  }, numeric(1L)))
  raw_observed_ss <- sum(vapply(seq_along(blocks), function(k) {
    sum(blocks[[k]][mask[[k]]]^2)
  }, numeric(1L)))

  expect_equal(
    diag$prediction_error,
    fit$missing$convergence$objective / observed_ss,
    tolerance = 1e-8
  )
  expect_gt(abs(diag$prediction_error -
                  fit$missing$convergence$objective / raw_observed_ss),
            1e-4)
})
