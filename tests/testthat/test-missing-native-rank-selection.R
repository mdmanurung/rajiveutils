test_that("native_cv rank selection runs only in native missing mode", {
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )

  expect_error(
    Rajive(blocks, c(2L, 2L), joint_rank = "native_cv"),
    class = "rajiveplus_invalid_input"
  )
})

test_that("native_cv records candidate diagnostics and selected rank", {
  set.seed(8801)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[2, 3] <- FALSE

  fit <- Rajive(blocks, c(2L, 2L), missing = "native",
                mask = mask, joint_rank = "native_cv",
                missing_control = rajive_missing_control(rank_candidates = 0:2))

  expect_s3_class(fit, "rajive_incomplete")
  expect_true(is.data.frame(fit$missing$rank_diagnostics))
  expect_true(fit$joint_rank %in% 0:2)
  expect_equal(fit$missing$selected_rank, fit$joint_rank)
})

test_that("native joint_rank NA defaults to native rank selection", {
  set.seed(8802)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block2[3, 2] <- FALSE

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask)

  expect_s3_class(fit, "rajive_incomplete")
  expect_true(is.data.frame(fit$missing$rank_diagnostics))
  expect_true(0L %in% fit$missing$rank_diagnostics$joint_rank)
  expect_equal(fit$missing$selected_rank, fit$joint_rank)
})

test_that("fixed native ranks do not run automatic diagnostics", {
  set.seed(8803)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[2, 3] <- FALSE

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L)

  expect_null(fit$missing$rank_diagnostics)
  expect_null(fit$missing$selected_rank)
  expect_equal(fit$joint_rank, 1L)
})

test_that("invalid native rank-selection strings are rejected", {
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )

  expect_error(
    Rajive(blocks, c(2L, 2L), missing = "native", joint_rank = "bad_rank"),
    class = "rajiveplus_invalid_input"
  )
})
