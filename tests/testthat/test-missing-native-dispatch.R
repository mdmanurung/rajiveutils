test_that("Rajive complete-data default remains a plain rajive fit", {
  set.seed(8401)
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )

  fit <- Rajive(blocks, c(2L, 2L), joint_rank = 1L,
                n_wedin_samples = NA, n_rand_dir_samples = NA)

  expect_s3_class(fit, "rajive")
  expect_false(inherits(fit, "rajive_incomplete"))
})

test_that("Rajive missing error mode rejects NAs", {
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  blocks$block1[2, 3] <- NA_real_

  expect_error(
    Rajive(blocks, c(2L, 2L), missing = "error"),
    class = "rajiveplus_invalid_input"
  )
})

test_that("Rajive native mode returns a rajive_incomplete fit", {
  set.seed(8402)
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  blocks$block1[2, 3] <- NA_real_

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", joint_rank = 1L)

  expect_s3_class(fit, "rajive_incomplete")
  expect_s3_class(fit, "rajive")
  expect_equal(fit$joint_rank, 1L)
})

test_that("Rajive native mode routes finite blocks when explicit mask is supplied", {
  set.seed(8403)
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block2[4, ] <- FALSE

  fit <- Rajive(blocks, c(2L, 2L), missing = "native",
                mask = mask, joint_rank = 1L)

  expect_s3_class(fit, "rajive_incomplete")
  expect_equal(nrow(get_joint_scores(fit)), 12L)
  expect_equal(get_joint_rank(fit), 1L)
})

test_that("fully observed native fits warn when preprocessing is requested", {
  set.seed(8404)
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))

  expect_warning(
    Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
           joint_rank = 1L,
           missing_control = rajive_missing_control(center = TRUE)),
    class = "rajiveplus_preprocess_skipped"
  )
})
