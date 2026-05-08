# Tests for visualization-related remediation items (W-M6, W-M9).

test_that("assess_stability components is permutation/sign invariant via Procrustes", {
  set.seed(777)
  n <- 30

  # Build reference scores with two components.
  ref_scores <- cbind(
    seq_len(n),
    rev(seq_len(n))
  )
  ref_scores <- scale(ref_scores)

  ajive_output <- structure(
    list(joint_scores = ref_scores, joint_rank = 2),
    class = "rajive"
  )

  # Encode row IDs in first column so mocked Rajive can reconstruct sampled idx.
  X1 <- cbind(seq_len(n), matrix(rnorm(n * 3), n, 3))
  X2 <- cbind(seq_len(n), matrix(rnorm(n * 2), n, 2))
  blocks <- list(X1, X2)

  # Mock Rajive: returns swapped + sign-flipped components for bootstrap sample.
  fake_rajive <- function(b_list, initial_signal_ranks) {
    idx <- as.integer(round(b_list[[1]][, 1]))
    bs <- ref_scores[idx, c(2, 1), drop = FALSE]
    bs[, 1] <- -bs[, 1]
    list(joint_scores = bs, joint_rank = 2)
  }

  out <- testthat::with_mocked_bindings(
    assess_stability(
      ajive_output = ajive_output,
      blocks = blocks,
      initial_signal_ranks = c(2, 2),
      target = "components",
      B = 10,
      sample_frac = 0.75
    ),
    Rajive = fake_rajive
  )

  expect_true(all(out$mean_correlation > 0.95))
})

test_that("associate_components(survival) rejects invalid status type", {
  skip_if_not_installed("survival")

  n <- 20
  ajive_output <- structure(
    list(joint_scores = matrix(rnorm(n * 2), n, 2), joint_rank = 2),
    class = "rajive"
  )

  metadata <- data.frame(
    time = rexp(n, rate = 0.1),
    status = sample(c("alive", "dead"), n, replace = TRUE)
  )

  expect_error(
    associate_components(
      ajive_output = ajive_output,
      metadata = metadata,
      mode = "survival",
      variable = "time",
      time_col = "time",
      status_col = "status"
    ),
    class = "rajiveplus_invalid_input"
  )
})
