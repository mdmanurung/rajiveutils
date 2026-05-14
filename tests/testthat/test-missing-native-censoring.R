test_that("native diagnostics state generic observed-entry assumptions", {
  set.seed(9001)
  blocks <- list(
    block1 = matrix(rnorm(50), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[1, 1] <- FALSE

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L)
  diag <- get_missing_diagnostics(fit)

  expect_match(diag$assumptions, "MCAR|MAR")
})

test_that("censoring metadata is accepted and reported", {
  set.seed(9002)
  blocks <- list(
    block1 = matrix(rnorm(50), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block2[3, 2] <- FALSE
  control <- rajive_missing_control(
    censoring = list(lower_limit = list(block2 = rep(-1, ncol(blocks$block2))))
  )

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L, missing_control = control)
  diag <- get_missing_diagnostics(fit)

  expect_true(diag$censoring$enabled)
  expect_equal(diag$censoring$blocks, "block2")
})

test_that("sensitivity metadata is recorded without fabricated deltas", {
  set.seed(9003)
  blocks <- list(
    block1 = matrix(rnorm(50), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[1, 1] <- FALSE
  control <- rajive_missing_control(sensitivity = c("mar", "censor_low"))

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L, missing_control = control)
  diag <- get_missing_diagnostics(fit)

  expect_equal(diag$sensitivity$assumptions, c("mar", "censor_low"))
  expect_false(diag$sensitivity$implemented)
  expect_match(diag$sensitivity$note, "not yet implemented")
})
