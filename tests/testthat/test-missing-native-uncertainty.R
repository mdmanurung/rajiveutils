test_that("bootstrap missing-data uncertainty stores aligned native refits", {
  set.seed(8901)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[1, 1] <- FALSE

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L, uncertainty = "bootstrap",
                missing_control = rajive_missing_control(n_refits = 3L))
  unc <- get_missing_uncertainty(fit)

  expect_equal(unc$method, "bootstrap")
  expect_equal(length(unc$refits), 3L)
  max_dev <- max(vapply(unc$refits, function(r) {
    max(abs(r$joint_scores[, 1] - fit$joint_scores[, 1]))
  }, numeric(1L)))
  expect_gt(max_dev, 1e-4)
})

test_that("mi missing-data uncertainty stores plausible completed analyses", {
  set.seed(8902)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block2[2, ] <- FALSE

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L, uncertainty = "mi",
                missing_control = rajive_missing_control(n_refits = 2L))
  unc <- get_missing_uncertainty(fit)

  expect_equal(unc$method, "mi")
  expect_equal(length(unc$refits), 2L)
  expect_true(is.data.frame(unc$score_summary))
  expect_match(unc$adapter_note, "rajive_ci")
})

test_that("mi uncertainty refits do not warn about skipped native preprocessing", {
  set.seed(8903)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[1, 1] <- FALSE
  control <- rajive_missing_control(center = TRUE, scale = TRUE,
                                    normalize = TRUE, n_refits = 2L)

  fit <- NULL
  expect_warning(
    fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                  joint_rank = 1L, uncertainty = "mi",
                  missing_control = control,
                  identifiability_norm = "l2"),
    regexp = NA
  )

  expect_equal(get_missing_uncertainty(fit)$method, "mi")
})
