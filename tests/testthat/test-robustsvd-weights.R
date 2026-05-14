test_that("RobRSVD.all without weights preserves the unweighted result", {
  set.seed(8101)
  X <- matrix(rnorm(40), nrow = 8)

  expected <- rajiveplus:::RobRSVD.all(X, nrank = 3)
  got <- rajiveplus:::RobRSVD.all(X, nrank = 3, weights = NULL)

  expect_equal(got$d, expected$d, tolerance = 1e-10)
  expect_equal(got$u, expected$u, tolerance = 1e-10)
  expect_equal(got$v, expected$v, tolerance = 1e-10)
})

test_that("RobRSVD.all all-one weights match no weights", {
  set.seed(8102)
  X <- matrix(rnorm(54), nrow = 9)
  W <- matrix(1, nrow = nrow(X), ncol = ncol(X))

  expected <- rajiveplus:::RobRSVD.all(X, nrank = 2)
  got <- rajiveplus:::RobRSVD.all(X, nrank = 2, weights = W)

  expect_equal(got$d, expected$d, tolerance = 1e-10)
  expect_equal(got$u, expected$u, tolerance = 1e-10)
  expect_equal(got$v, expected$v, tolerance = 1e-10)
})

test_that("weighted robust SVD ignores corrupted masked cells", {
  set.seed(8103)
  X <- matrix(rnorm(60), nrow = 10)
  W <- matrix(1, nrow = nrow(X), ncol = ncol(X))
  W[3, 4] <- 0
  W[8, 2] <- 0

  X_corrupt <- X
  X_corrupt[3, 4] <- 1e9
  X_corrupt[8, 2] <- -1e9

  a <- rajiveplus:::RobRSVD.all(X, nrank = 2, weights = W)
  b <- rajiveplus:::RobRSVD.all(X_corrupt, nrank = 2, weights = W)

  expect_equal(a$d, b$d, tolerance = 1e-8)
  expect_equal(abs(crossprod(a$u, b$u)), diag(2), tolerance = 1e-6)
  expect_equal(abs(crossprod(a$v, b$v)), diag(2), tolerance = 1e-6)
})

test_that("weighted robust SVD handles fully masked rows and columns", {
  set.seed(8104)
  X <- matrix(rnorm(48), nrow = 8)
  W <- matrix(1, nrow = nrow(X), ncol = ncol(X))
  W[2, ] <- 0
  W[, 5] <- 0

  got <- rajiveplus:::RobRSVD.all(X, nrank = 2, weights = W)

  expect_equal(dim(got$u), c(8L, 2L))
  expect_equal(dim(got$v), c(6L, 2L))
  expect_true(all(is.finite(got$d)))
  expect_equal(got$u[2, ], c(0, 0), tolerance = 1e-12)
  expect_equal(got$v[5, ], c(0, 0), tolerance = 1e-12)
})

test_that("weighted robust SVD can shrink singular values", {
  set.seed(8105)
  X <- matrix(rnorm(54), nrow = 9)
  W <- matrix(1, nrow = nrow(X), ncol = ncol(X))

  base <- rajiveplus:::RobRSVD.all(X, nrank = 2, weights = W)
  shrunk <- rajiveplus:::RobRSVD.all(X, nrank = 2, weights = W,
                                     shrinkage = 0.25)

  expect_equal(shrunk$d, pmax(base$d - 0.25, 0), tolerance = 1e-10)
  expect_lt(norm(rajiveplus:::svd_reconstruction(shrunk), "F"),
            norm(rajiveplus:::svd_reconstruction(base), "F"))
})

test_that("missMDA-style shrinkage uses discarded singular values", {
  sv <- c(5, 3, 2, 1)
  got <- rajiveplus:::.shrink_singular_values(sv, rank = 2,
                                              shrinkage = "missmda",
                                              shrinkage_coeff = 1,
                                              n_rows = 10,
                                              n_cols = 4)
  denom <- (10 - 1) * 4 - (10 - 1) * 2 - 4 * 2 + 2^2
  sigma2 <- 10 * 4 / min(4, 10 - 1) * sum(c(2, 1)^2 / denom)
  sigma2 <- min(sigma2, 2^2)
  expect_equal(got, c((5^2 - sigma2) / 5, (3^2 - sigma2) / 3))
})

test_that("RobRSVD.all validates weight dimensions", {
  X <- matrix(rnorm(20), nrow = 5)
  W <- matrix(1, nrow = 4, ncol = 4)

  expect_error(
    rajiveplus:::RobRSVD.all(X, nrank = 2, weights = W),
    class = "rajiveplus_invalid_input"
  )
})
