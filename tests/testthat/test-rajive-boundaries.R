library(testthat)
library(rajiveplus)

# ---------------------------------------------------------------------------
# W-H1: get_sv_threshold() NA at boundary rank == length(d)
# ---------------------------------------------------------------------------

test_that("get_sv_threshold returns finite value at rank == length(d)", {
  # Red test from PLANS.md W-H1: was NA before fix because sv[rank+1] is NA.
  expect_false(is.na(rajiveplus:::get_sv_threshold(c(3, 2, 1), rank = 3)))
  # Concrete value: midpoint between d[3]=1 and implicit floor 0 = 0.5.
  expect_equal(rajiveplus:::get_sv_threshold(c(3, 2, 1), rank = 3), 0.5)
})

test_that("get_sv_threshold returns finite value at rank == length(d) length-1", {
  expect_false(is.na(rajiveplus:::get_sv_threshold(c(5), rank = 1)))
  expect_equal(rajiveplus:::get_sv_threshold(c(5), rank = 1), 2.5)
})

test_that("get_sv_threshold matches midpoint formula in interior", {
  # Regression: interior entries must be unchanged.
  d <- c(10, 5, 2, 1)
  expect_equal(rajiveplus:::get_sv_threshold(d, rank = 1), 0.5 * (10 + 5))
  expect_equal(rajiveplus:::get_sv_threshold(d, rank = 2), 0.5 * (5  + 2))
  expect_equal(rajiveplus:::get_sv_threshold(d, rank = 3), 0.5 * (2  + 1))
})

test_that("Rajive does not error when initial_signal_ranks equals min(dim(block))", {
  # Trigger the boundary path: 5x4 blocks, initial_signal_ranks = c(4, 4),
  # so svd_ranks = pmin(4+1, 4) = c(4, 4) and get_sv_threshold is called at
  # rank == length(d).
  set.seed(1)
  b1 <- matrix(rnorm(5 * 4), 5, 4)
  b2 <- matrix(rnorm(5 * 4), 5, 4)
  fit <- expect_warning(
    Rajive(list(b1, b2), initial_signal_ranks = c(4, 4), joint_rank = 1),
    class = "rajiveplus_underdetermined"
  )
  expect_s3_class(fit, "rajive")
})

# ---------------------------------------------------------------------------
# W-M3: truncate_svd(rank = 0) returns zero-column matrices
# ---------------------------------------------------------------------------

test_that("truncate_svd(rank = 0) returns zero-column matrices", {
  decomp <- list(
    u = matrix(rnorm(10), 5, 2),
    d = c(2, 1),
    v = matrix(rnorm(8), 4, 2)
  )
  out <- rajiveplus:::truncate_svd(decomp, rank = 0)
  expect_equal(ncol(out$u), 0L)
  expect_equal(ncol(out$v), 0L)
  expect_equal(length(out$d), 0L)
  expect_equal(nrow(out$u), 5L)
  expect_equal(nrow(out$v), 4L)
})

test_that("truncate_svd(rank = 0) reconstructs to zero matrix", {
  decomp <- list(
    u = matrix(rnorm(10), 5, 2),
    d = c(2, 1),
    v = matrix(rnorm(6), 3, 2)
  )
  out   <- rajiveplus:::truncate_svd(decomp, rank = 0)
  recon <- rajiveplus:::svd_reconstruction(out)
  expect_equal(dim(recon), c(5L, 3L))
  expect_true(all(recon == 0))
})

test_that("truncate_svd positive rank still works after W-M3 fix", {
  decomp <- list(
    u = matrix(rnorm(15), 5, 3),
    d = c(5, 3, 1),
    v = matrix(rnorm(9), 3, 3)
  )
  out <- rajiveplus:::truncate_svd(decomp, rank = 2)
  expect_equal(ncol(out$u), 2L)
  expect_equal(ncol(out$v), 2L)
  expect_equal(length(out$d), 2L)
})

# ---------------------------------------------------------------------------
# Regression: zero-rank path in get_joint_decomposition_robustH still correct
# ---------------------------------------------------------------------------

test_that("get_joint_decomposition_robustH zero-rank still produces zero-column u", {
  set.seed(42)
  X <- matrix(rnorm(20 * 10), 20, 10)
  # Pass null joint_scores (0-column matrix) to force zero-rank path.
  jd <- rajiveplus:::get_joint_decomposition_robustH(
    X, joint_scores = matrix(0, 20, 0), full = FALSE
  )
  expect_equal(ncol(jd$u), 0L)
  expect_equal(ncol(jd$v), 0L)
  expect_equal(length(jd$d), 0L)
})
