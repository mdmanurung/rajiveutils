# Regression tests for robust deflation behavior (W-H3/W-M7).

reconstruct_svd <- function(svd_list, r) {
  U <- svd_list$u[, seq_len(r), drop = FALSE]
  d <- svd_list$d[seq_len(r)]
  V <- svd_list$v[, seq_len(r), drop = FALSE]
  U %*% (diag(d, nrow = r, ncol = r) %*% t(V))
}

test_that("RobRSVD.all rank-2 fit improves over rank-1 on rank-2 signal", {
  skip_on_cran()

  set.seed(202)
  n <- 40
  p <- 30

  U0 <- qr.Q(qr(matrix(rnorm(n * 2), n, 2)))
  V0 <- qr.Q(qr(matrix(rnorm(p * 2), p, 2)))
  D0 <- diag(c(7, 4), 2, 2)
  X <- U0 %*% D0 %*% t(V0) + matrix(rnorm(n * p, sd = 0.05), n, p)

  # Inject sparse outliers to force robust iterations to matter.
  out_idx <- sample.int(n * p, 30)
  X[out_idx] <- X[out_idx] + rnorm(length(out_idx), sd = 8)

  fit2 <- rajiveplus:::RobRSVD.all(X, nrank = 2)
  fit1 <- rajiveplus:::RobRSVD.all(X, nrank = 1)

  Xhat2 <- reconstruct_svd(fit2, 2)
  Xhat1 <- reconstruct_svd(fit1, 1)

  err2 <- norm(X - Xhat2, type = "F")
  err1 <- norm(X - Xhat1, type = "F")

  expect_true(all(is.finite(fit2$d[1:2])))
  expect_gt(fit2$d[2], 0)
  expect_lt(err2, err1)
})

test_that("RobRSVD.all with explicit svdinit equals lazy svdinit pathway", {
  set.seed(303)
  X <- matrix(rnorm(25 * 20), 25, 20)

  fit_lazy <- rajiveplus:::RobRSVD.all(X, nrank = 2)
  fit_exp  <- rajiveplus:::RobRSVD.all(X, nrank = 2, svdinit = svd(X))

  expect_equal(fit_lazy$d[1:2], fit_exp$d[1:2], tolerance = 1e-8)
})
