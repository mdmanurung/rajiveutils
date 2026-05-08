# Regression tests for variance explained computation (W-M2).

make_mock_decomp <- function(blocks) {
  # Build 3*K decomposition list in AJIVE order: I_k, J_k, E_k.
  # Put `d` last in each list to guard against positional [[1]] bugs.
  list(
    list(u = matrix(10, nrow(blocks[[1]]), 2), v = matrix(0, ncol(blocks[[1]]), 2), d = c(0.8, 0.5)),
    list(u = matrix(99, nrow(blocks[[1]]), 2), v = matrix(0, ncol(blocks[[1]]), 2), d = c(1.1, 0.7)),
    list(u = matrix(0, nrow(blocks[[1]]), 1),  v = matrix(0, ncol(blocks[[1]]), 1), d = 0),
    list(u = matrix(10, nrow(blocks[[2]]), 2), v = matrix(0, ncol(blocks[[2]]), 2), d = c(0.9, 0.4)),
    list(u = matrix(99, nrow(blocks[[2]]), 2), v = matrix(0, ncol(blocks[[2]]), 2), d = c(1.0, 0.6)),
    list(u = matrix(0, nrow(blocks[[2]]), 1),  v = matrix(0, ncol(blocks[[2]]), 1), d = 0)
  )
}

test_that("showVarExplained_robust uses named 'd' singular values", {
  set.seed(100)
  blocks <- list(
    matrix(rnorm(5 * 4), 5, 4),
    matrix(rnorm(5 * 3), 5, 3)
  )

  fit <- list(block_decomps = make_mock_decomp(blocks))
  got <- showVarExplained_robust(fit, blocks)

  # If positional [[1]] were used here, huge 'u' entries would inflate > 1.
  expect_true(all(got$Joint <= 1 + 1e-12))
  expect_true(all(got$Indiv <= 1 + 1e-12))
})

test_that("showVarExplained_robust matches sum(d^2)/||X||_F^2", {
  set.seed(101)
  blocks <- list(
    matrix(rnorm(6 * 5), 6, 5),
    matrix(rnorm(6 * 4), 6, 4)
  )

  fit <- list(block_decomps = make_mock_decomp(blocks))
  got <- showVarExplained_robust(fit, blocks)

  expected_joint <- c(
    sum(c(1.1, 0.7)^2) / norm(blocks[[1]], "F")^2,
    sum(c(1.0, 0.6)^2) / norm(blocks[[2]], "F")^2
  )
  expected_indiv <- c(
    sum(c(0.8, 0.5)^2) / norm(blocks[[1]], "F")^2,
    sum(c(0.9, 0.4)^2) / norm(blocks[[2]], "F")^2
  )

  expect_equal(got$Joint, expected_joint)
  expect_equal(got$Indiv, expected_indiv)
})

test_that("showVarExplained_robust components sum to 1 per block", {
  set.seed(102)
  blocks <- list(
    matrix(rnorm(8 * 5), 8, 5),
    matrix(rnorm(8 * 4), 8, 4)
  )

  fit <- list(block_decomps = make_mock_decomp(blocks))
  got <- showVarExplained_robust(fit, blocks)

  expect_equal(got$Joint + got$Indiv + got$Resid,
               rep(1, length(blocks)),
               tolerance = 1e-10)
})
