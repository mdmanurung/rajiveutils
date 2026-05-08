library(testthat)
library(rajiveplus)

# ---------------------------------------------------------------------------
# Reference implementations (pure R, verbatim copy of original RobustSVD.R)
# Used to verify that Rcpp versions produce identical numerical results.
# ---------------------------------------------------------------------------

.RobRSVD1_R <- function(data, huberk = 1.345, niter = 1000,
                        tol = 1e-05, sinit, uinit, vinit) {
  size_data <- c(dim(data))
  m <- size_data[1]
  n <- size_data[2]
  sold <- sinit
  vold <- vinit
  uold <- sold * uinit
  Appold <- uold %*% t(vold)
  Rmat <- data - Appold
  Rvec <- c(Rmat)
  mysigma <- median(abs(Rvec)) / 0.675
  if (!is.finite(mysigma) || mysigma <= 0) {
    mysigma <- sqrt(.Machine$double.eps)
  }
  iter <- 1
  localdiff <- 9999
  while (localdiff > tol & iter < niter) {
    # W-M7 reference parity: robust scale is updated every iteration.
    Rvec <- c(Rmat)
    mysigma <- median(abs(Rvec)) / 0.675
    if (!is.finite(mysigma) || mysigma <= 0) {
      mysigma <- sqrt(.Machine$double.eps)
    }

    Wmat <- huberk / abs(Rmat / mysigma)
    Wmat[Wmat > 1] <- 1
    uterm1 <- diag(colSums(diag(c(vold^2)) %*% t(Wmat))) +
      (2 * mysigma^2) * (c(t(vold) %*% (diag(n)) %*% vold) * (diag(m)) - diag(sum(vold^2), m))
    uterm2 <- (Wmat * data) %*% vold
    unew <- solve(uterm1) %*% uterm2
    vterm1 <- diag(colSums(diag(c(unew^2)) %*% Wmat)) +
      (2 * mysigma^2) * (c(t(unew) %*% (diag(m)) %*% unew) * (diag(n)) - diag(sum(unew^2), n))
    vterm2 <- t(Wmat * data) %*% unew
    vnew <- solve(vterm1) %*% vterm2
    Appnew <- unew %*% t(vnew)
    Rmat <- data - Appnew
    localdiff <- max(abs(Appnew - Appold))
    Appold <- Appnew
    uold <- sqrt(sum(vnew^2)) * unew
    vold <- vnew / sqrt(sum(vnew^2))
    iter <- iter + 1
  }
  v <- vold
  s <- sqrt(sum(uold^2))
  u <- uold / sqrt(sum(uold^2))
  list(s = s, u = u, v = v)
}

.RobRSVD_all_R <- function(data, nrank = min(dim(data)), svdinit = svd(data)) {
  data.svd1 <- .RobRSVD1_R(data, sinit = svdinit$d[1],
                            uinit = svdinit$u[, 1], vinit = svdinit$v[, 1])
  d   <- data.svd1$s
  u   <- data.svd1$u
  v   <- data.svd1$v
  Red <- d * u %*% t(v)
  Rm  <- min(min(dim(data)), nrank)
  for (i in 1:(Rm - 1)) {
    # W-H3 reference parity: warm-start from leading SVD triplet
    # of the current residual, not from original svdinit.
    residual <- data - Red
    svd_res  <- svd(residual)
    data.svd1 <- .RobRSVD1_R(residual,
                             sinit = svd_res$d[1],
                             uinit = svd_res$u[, 1],
                             vinit = svd_res$v[, 1])
    d   <- c(d, data.svd1$s)
    u   <- cbind(u, data.svd1$u)
    v   <- cbind(v, data.svd1$v)
    Red <- u %*% diag(d) %*% t(v)
  }
  list(d = d, u = u, v = v)
}

# Helper: compare rank-r reconstruction regardless of column sign flips.
# sign-invariant: project each singular vector to a rank-1 outer product and
# compare the outer products (u_j u_j^T and v_j v_j^T).
expect_svd_equal <- function(r_res, cpp_res, nrank, tol = 1e-7) {
  # singular values — always non-negative, no sign issue
  expect_equal(as.numeric(r_res$d), as.numeric(cpp_res$d),
               tolerance = tol,
               label = "singular values")

  for (j in seq_len(nrank)) {
    uj_r   <- as.numeric(r_res$u[, j])
    uj_cpp <- as.numeric(cpp_res$u[, j])
    vj_r   <- as.numeric(r_res$v[, j])
    vj_cpp <- as.numeric(cpp_res$v[, j])

    # Projection matrices are sign-invariant: u u^T
    expect_equal(tcrossprod(matrix(uj_r)),
                 tcrossprod(matrix(uj_cpp)),
                 tolerance = tol,
                 label = paste0("u projection, component ", j))
    expect_equal(tcrossprod(matrix(vj_r)),
                 tcrossprod(matrix(vj_cpp)),
                 tolerance = tol,
                 label = paste0("v projection, component ", j))

    # Full rank-1 reconstruction: d * u * v^T (sign in u and v must be consistent)
    expect_equal(r_res$d[j]   * outer(uj_r, vj_r),
                 cpp_res$d[j] * outer(uj_cpp, vj_cpp),
                 tolerance = tol,
                 label = paste0("rank-1 reconstruction, component ", j))
  }
}


# ===========================================================================
# 1. RobRSVD1_cpp vs pure-R: rank-1 component only
# ===========================================================================

test_that("RobRSVD1_cpp matches R: tall matrix (rows > cols)", {
  set.seed(101)
  m <- 30; n <- 10
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  r_res  <- .RobRSVD1_R(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])
  cpp_res <- RobRSVD1_cpp(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])

  expect_equal(r_res$s, cpp_res$s, tolerance = 1e-7)
  # sign-invariant projection matrices
  expect_equal(tcrossprod(matrix(r_res$u)),   tcrossprod(matrix(cpp_res$u)),   tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$v)),   tcrossprod(matrix(cpp_res$v)),   tolerance = 1e-7)
  # sign-consistent reconstruction
  expect_equal(r_res$s   * outer(as.numeric(r_res$u),   as.numeric(r_res$v)),
               cpp_res$s * outer(as.numeric(cpp_res$u), as.numeric(cpp_res$v)), tolerance = 1e-7)
})

test_that("RobRSVD1_cpp matches R: wide matrix (cols > rows)", {
  set.seed(102)
  m <- 10; n <- 25
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  r_res  <- .RobRSVD1_R(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])
  cpp_res <- RobRSVD1_cpp(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])

  expect_equal(r_res$s, cpp_res$s, tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$u)),   tcrossprod(matrix(cpp_res$u)),   tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$v)),   tcrossprod(matrix(cpp_res$v)),   tolerance = 1e-7)
  expect_equal(r_res$s   * outer(as.numeric(r_res$u),   as.numeric(r_res$v)),
               cpp_res$s * outer(as.numeric(cpp_res$u), as.numeric(cpp_res$v)), tolerance = 1e-7)
})

test_that("RobRSVD1_cpp matches R: square matrix", {
  set.seed(103)
  m <- 15
  X      <- matrix(rnorm(m * m), m, m)
  sv0    <- svd(X)
  r_res  <- .RobRSVD1_R(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])
  cpp_res <- RobRSVD1_cpp(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])

  expect_equal(r_res$s, cpp_res$s, tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$u)),   tcrossprod(matrix(cpp_res$u)),   tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$v)),   tcrossprod(matrix(cpp_res$v)),   tolerance = 1e-7)
  expect_equal(r_res$s   * outer(as.numeric(r_res$u),   as.numeric(r_res$v)),
               cpp_res$s * outer(as.numeric(cpp_res$u), as.numeric(cpp_res$v)), tolerance = 1e-7)
})

test_that("RobRSVD1_cpp matches R: near-rank-1 matrix with small noise", {
  set.seed(104)
  m <- 20; n <- 15
  # Low-rank signal + small noise
  u_true <- rnorm(m); u_true <- u_true / sqrt(sum(u_true^2))
  v_true <- rnorm(n); v_true <- v_true / sqrt(sum(v_true^2))
  X      <- 10 * outer(u_true, v_true) + 0.01 * matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  r_res  <- .RobRSVD1_R(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])
  cpp_res <- RobRSVD1_cpp(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])

  expect_equal(r_res$s, cpp_res$s, tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$u)),   tcrossprod(matrix(cpp_res$u)),   tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$v)),   tcrossprod(matrix(cpp_res$v)),   tolerance = 1e-7)
})

test_that("RobRSVD1_cpp matches R: custom huberk, niter, tol", {
  set.seed(105)
  m <- 18; n <- 12
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  r_res   <- .RobRSVD1_R(X, huberk = 1.0, niter = 200, tol = 1e-4,
                          sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])
  cpp_res <- RobRSVD1_cpp(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1],
                          huberk = 1.0, niter = 200L, tol = 1e-4)

  expect_equal(r_res$s, cpp_res$s, tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$u)),   tcrossprod(matrix(cpp_res$u)),   tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$v)),   tcrossprod(matrix(cpp_res$v)),   tolerance = 1e-7)
})

test_that("RobRSVD1_cpp matches R: matrix with outliers", {
  set.seed(106)
  m <- 25; n <- 10
  X    <- matrix(rnorm(m * n), m, n)
  # Inject outliers in a few cells
  X[sample(m * n, size = 10)] <- rnorm(10, mean = 0, sd = 20)
  sv0    <- svd(X)
  r_res  <- .RobRSVD1_R(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])
  cpp_res <- RobRSVD1_cpp(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])

  expect_equal(r_res$s, cpp_res$s, tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$u)),   tcrossprod(matrix(cpp_res$u)),   tolerance = 1e-7)
  expect_equal(tcrossprod(matrix(r_res$v)),   tcrossprod(matrix(cpp_res$v)),   tolerance = 1e-7)
})


# ===========================================================================
# 2. RobRSVD_all_cpp vs pure-R: multi-rank
# ===========================================================================

test_that("RobRSVD_all_cpp matches R: rank-3 tall matrix", {
  set.seed(201)
  m <- 40; n <- 20; nrank <- 3
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  r_res   <- .RobRSVD_all_R(X, nrank = nrank, svdinit = sv0)
  cpp_res <- RobRSVD_all_cpp(X, nrank = nrank,
                              sinit1 = sv0$d[1], uinit1 = sv0$u[, 1], vinit1 = sv0$v[, 1])

  expect_svd_equal(r_res, cpp_res, nrank = nrank, tol = 1e-7)
})

test_that("RobRSVD_all_cpp matches R: rank-2 wide matrix", {
  set.seed(202)
  m <- 15; n <- 30; nrank <- 2
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  r_res   <- .RobRSVD_all_R(X, nrank = nrank, svdinit = sv0)
  cpp_res <- RobRSVD_all_cpp(X, nrank = nrank,
                              sinit1 = sv0$d[1], uinit1 = sv0$u[, 1], vinit1 = sv0$v[, 1])

  expect_svd_equal(r_res, cpp_res, nrank = nrank, tol = 1e-7)
})

test_that("RobRSVD_all_cpp matches R: rank-4 with outliers", {
  set.seed(203)
  m <- 50; n <- 25; nrank <- 4
  X    <- matrix(rnorm(m * n), m, n)
  X[sample(m * n, size = 20)] <- rnorm(20, sd = 15)
  sv0    <- svd(X)
  r_res   <- .RobRSVD_all_R(X, nrank = nrank, svdinit = sv0)
  cpp_res <- RobRSVD_all_cpp(X, nrank = nrank,
                              sinit1 = sv0$d[1], uinit1 = sv0$u[, 1], vinit1 = sv0$v[, 1])

  expect_svd_equal(r_res, cpp_res, nrank = nrank, tol = 1e-7)
})

test_that("RobRSVD_all_cpp matches R: rank-1 (single component)", {
  # Note: R's `for(i in 1:(Rm-1))` when Rm=1 evaluates `1:0 = c(1,0)` and
  # iterates *twice*, producing 3 components — an unintentional R idiom.
  # The C++ version correctly returns 1 component for nrank=1.
  # We therefore only verify the structural correctness of the C++ result here
  # and do NOT compare against the R reference for this edge case.
  # In the RaJIVE pipeline nrank >= 2 always holds in practice.
  set.seed(204)
  m <- 20; n <- 12; nrank <- 1
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  cpp_res <- RobRSVD_all_cpp(X, nrank = nrank,
                              sinit1 = sv0$d[1], uinit1 = sv0$u[, 1], vinit1 = sv0$v[, 1])

  expect_length(cpp_res$d, 1L)
  expect_equal(dim(cpp_res$u), c(m, 1L))
  expect_equal(dim(cpp_res$v), c(n, 1L))
  expect_gt(cpp_res$d[1], 0)
})


# ===========================================================================
# 3. get_svd_robustH end-to-end (the R wrapper used throughout the pipeline)
# ===========================================================================

test_that("get_svd_robustH matches pure-R RobRSVD.all: no rank arg", {
  set.seed(301)
  m <- 25; n <- 10
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  r_res   <- .RobRSVD_all_R(X, nrank = min(dim(X)), svdinit = sv0)
  cpp_res <- rajiveplus:::get_svd_robustH(X)

  expect_svd_equal(r_res, cpp_res, nrank = n, tol = 1e-7)
})

test_that("get_svd_robustH matches pure-R RobRSVD.all: with rank arg", {
  set.seed(302)
  m <- 30; n <- 20; nrank <- 5
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  r_res   <- .RobRSVD_all_R(X, nrank = nrank, svdinit = sv0)
  cpp_res <- rajiveplus:::get_svd_robustH(X, rank = nrank)

  expect_svd_equal(r_res, cpp_res, nrank = nrank, tol = 1e-7)
})


# ===========================================================================
# 4. Return structure validation
# ===========================================================================

test_that("RobRSVD1_cpp returns list with s, u, v of correct dimensions", {
  set.seed(401)
  m <- 12; n <- 8
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  res    <- RobRSVD1_cpp(X, sinit = sv0$d[1], uinit = sv0$u[, 1], vinit = sv0$v[, 1])

  expect_named(res, c("s", "u", "v"))
  expect_length(res$s, 1L)
  expect_length(res$u, m)
  expect_length(res$v, n)
  expect_true(is.numeric(res$s))
  expect_gt(res$s, 0)
})

test_that("RobRSVD_all_cpp returns list with d, u, v of correct dimensions", {
  set.seed(402)
  m <- 20; n <- 10; nrank <- 3
  X      <- matrix(rnorm(m * n), m, n)
  sv0    <- svd(X)
  res    <- RobRSVD_all_cpp(X, nrank = nrank,
                             sinit1 = sv0$d[1], uinit1 = sv0$u[, 1], vinit1 = sv0$v[, 1])

  expect_named(res, c("d", "u", "v"))
  expect_length(res$d, nrank)
  expect_equal(dim(res$u), c(m, nrank))
  expect_equal(dim(res$v), c(n, nrank))
  expect_true(all(res$d > 0))
})
