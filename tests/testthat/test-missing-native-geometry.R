test_that("masked Frobenius norm ignores corrupted masked values", {
  X <- matrix(1:9, nrow = 3)
  W <- matrix(TRUE, nrow = 3, ncol = 3)
  W[2, 2] <- FALSE
  X_corrupt <- X
  X_corrupt[2, 2] <- 1e8

  expect_equal(
    rajiveplus:::.masked_frob_norm(X, W),
    rajiveplus:::.masked_frob_norm(X_corrupt, W),
    tolerance = 1e-12
  )
})

test_that("masked residual sum of squares equals observed-cell hand calculation", {
  X <- matrix(c(1, 2, 3, 4), nrow = 2)
  fitted <- matrix(c(1, 1, 1, 1), nrow = 2)
  W <- matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2)

  expect_equal(
    rajiveplus:::.masked_residual_ss(X, fitted, W),
    sum((X[W] - fitted[W])^2)
  )
})

test_that("masked variance explained partitions complete-data variance", {
  X <- matrix(c(2, 0, 0, 3), nrow = 2)
  J <- matrix(c(1, 0, 0, 1), nrow = 2)
  I <- matrix(c(1, 0, 0, 1), nrow = 2)
  W <- matrix(TRUE, nrow = 2, ncol = 2)

  got <- rajiveplus:::.masked_variance_explained(X, J, I, W)
  expect_equal(sum(got[c("Joint", "Indiv", "Resid")]), 1, tolerance = 1e-12)
  expect_equal(got[["observed_fraction"]], 1)
})

test_that("masked orthogonality detects non-orthogonal components", {
  W <- matrix(TRUE, nrow = 3, ncol = 2)
  A <- matrix(c(1, 0, 0, 1, 0, 0), nrow = 3)
  B <- A

  got <- rajiveplus:::.masked_orthogonality(A, B, W)
  expect_gt(got[["max_abs_crossprod"]], 0.9)
})
