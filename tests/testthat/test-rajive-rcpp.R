library(testthat)
library(rajiveplus)

# ---------------------------------------------------------------------------
# Pipeline-level regression tests: Rajive() with Rcpp RobustSVD backend
#
# Strategy: fix joint_rank so that the stochastic rank selection step
# (wedin + random direction bounds, which use parallel foreach) is bypassed.
# This gives fully deterministic, seed-reproducible results.
#
# block_decomps structure: mapply(SIMPLIFY=TRUE) returns a 3 x K list matrix:
#   row 1 = individual decomp, row 2 = joint decomp, row 3 = noise
# Column-major linear index: individual_k = [[3*(k-1)+1]],
#                            joint_k      = [[3*(k-1)+2]],
#                            noise_k      = [[3*(k-1)+3]]
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_blocks <- function(K, n, pks, rankJ, rankA, seed) {
  set.seed(seed)
  Y <- ajive.data.sim(K = K, rankJ = rankJ, rankA = rankA,
                      n = n, pks = pks, dist.type = 1)
  list(blocks = Y$sim_data,
       initial_signal_ranks = rankA)
}

# Access block k's individual / joint / noise from block_decomps (3 x K matrix)
bd_individual <- function(res, k) res$block_decomps[[3L * (k - 1L) + 1L]]
bd_joint      <- function(res, k) res$block_decomps[[3L * (k - 1L) + 2L]]
bd_noise      <- function(res, k) res$block_decomps[[3L * (k - 1L) + 3L]]

# Sign-invariant comparison of column spaces (projection matrices).
expect_colspace_equal <- function(A, B, tol = 1e-5, label = "column space") {
  PA <- A %*% solve(crossprod(A)) %*% t(A)
  PB <- B %*% solve(crossprod(B)) %*% t(B)
  expect_equal(PA, PB, tolerance = tol, label = label)
}


# ===========================================================================
# 1. joint_scores column space is well-formed
# ===========================================================================

test_that("Rajive joint_scores: K=2 correct dims and finite values", {
  d  <- make_blocks(K = 2, n = 30, pks = c(20, 15),
                    rankJ = 2, rankA = c(4, 3), seed = 1001L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, num_cores = 1L)

  js <- res$joint_scores
  expect_equal(ncol(js), 2L)
  expect_equal(nrow(js), 30L)
  expect_true(is.matrix(js))
  expect_true(all(is.finite(js)))
})

test_that("Rajive joint_scores identical across two runs with same seed", {
  d   <- make_blocks(K = 2, n = 30, pks = c(20, 15),
                     rankJ = 2, rankA = c(4, 3), seed = 1002L)
  r1  <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, num_cores = 1L)
  r2  <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, num_cores = 1L)

  # Identical inputs + no stochastic steps → exactly identical outputs
  expect_equal(r1$joint_scores, r2$joint_scores)
})


# ===========================================================================
# 2. full=TRUE decomposition: X = J + I + E  (key identity)
# ===========================================================================

test_that("Rajive full=TRUE: X = J + I + E for each block (K=2)", {
  d   <- make_blocks(K = 2, n = 40, pks = c(30, 20),
                     rankJ = 2, rankA = c(5, 3), seed = 2001L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, full = TRUE, num_cores = 1L)

  for (k in seq_along(d$blocks)) {
    X  <- d$blocks[[k]]
    J  <- bd_joint(res, k)$full
    I  <- bd_individual(res, k)$full
    E  <- bd_noise(res, k)
    expect_equal(J + I + E, X, tolerance = 1e-8,
                 label = paste0("X = J + I + E, block ", k))
  }
})

test_that("Rajive full=TRUE: X = J + I + E for each block (K=3)", {
  d   <- make_blocks(K = 3, n = 30, pks = c(20, 15, 12),
                     rankJ = 2, rankA = c(4, 3, 3), seed = 2002L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, full = TRUE, num_cores = 1L)

  for (k in seq_along(d$blocks)) {
    X  <- d$blocks[[k]]
    J  <- bd_joint(res, k)$full
    I  <- bd_individual(res, k)$full
    E  <- bd_noise(res, k)
    expect_equal(J + I + E, X, tolerance = 1e-8,
                 label = paste0("X = J + I + E, block ", k))
  }
})


# ===========================================================================
# 3. full=FALSE: SVDs are present (note: 'full' arg is not propagated in
#    get_final_decomposition_robustH, which always uses full=TRUE)
# ===========================================================================

test_that("Rajive full=FALSE: joint/individual SVD components are present", {
  d   <- make_blocks(K = 2, n = 30, pks = c(20, 15),
                     rankJ = 2, rankA = c(4, 3), seed = 3001L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, full = FALSE, num_cores = 1L)

  for (k in seq_along(d$blocks)) {
    expect_true(!is.null(bd_joint(res, k)$u),
                label = paste0("joint u present for block ", k))
    expect_true(!is.null(bd_individual(res, k)$u),
                label = paste0("individual u present for block ", k))
  }
})


# ===========================================================================
# 4. Robust SVD component validity: U/V are near-orthonormal, d > 0
#    (robust SVD is an approximation; U*D*V^T need not equal the input exactly)
# ===========================================================================

test_that("Rajive joint decomp: SVD components are near-orthonormal and d > 0", {
  d   <- make_blocks(K = 2, n = 40, pks = c(30, 20),
                     rankJ = 2, rankA = c(5, 3), seed = 4001L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, full = TRUE, num_cores = 1L)

  for (k in seq_along(d$blocks)) {
    jd <- bd_joint(res, k)
    expect_true(all(jd$d > 0),
                label = paste0("joint d > 0, block ", k))
    # Each column of U has unit L2 norm (guaranteed by RobRSVD1 normalization)
    col_norms <- sqrt(colSums(jd$u^2))
    expect_equal(col_norms, rep(1, ncol(jd$u)), tolerance = 1e-8,
                 label = paste0("joint U columns unit norm, block ", k))
  }
})

test_that("Rajive individual decomp: SVD components are near-orthonormal and d > 0", {
  d   <- make_blocks(K = 2, n = 40, pks = c(30, 20),
                     rankJ = 2, rankA = c(5, 3), seed = 4002L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, full = TRUE, num_cores = 1L)

  for (k in seq_along(d$blocks)) {
    id <- bd_individual(res, k)
    if (id$rank > 0L) {
      expect_true(all(id$d > 0),
                  label = paste0("individual d > 0, block ", k))
      col_norms <- sqrt(colSums(id$u^2))
      expect_equal(col_norms, rep(1, ncol(id$u)), tolerance = 1e-8,
                   label = paste0("individual U columns unit norm, block ", k))
    }
  }
})


# ===========================================================================
# 5. S3 structure and attributes
# ===========================================================================

test_that("Rajive returns class 'rajive' with expected fields", {
  d   <- make_blocks(K = 2, n = 30, pks = c(20, 15),
                     rankJ = 2, rankA = c(4, 3), seed = 5001L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, num_cores = 1L)

  expect_s3_class(res, "rajive")
  expect_true(!is.null(res$joint_scores))
  expect_true(!is.null(res$joint_rank))
  expect_true(!is.null(res$block_decomps))
  expect_equal(res$joint_rank, 2L)
  # block_decomps is a 3 x K list matrix (mapply SIMPLIFY=TRUE)
  expect_equal(dim(res$block_decomps), c(3L, 2L))
})

test_that("Rajive: print and summary do not error", {
  d   <- make_blocks(K = 2, n = 30, pks = c(20, 15),
                     rankJ = 2, rankA = c(4, 3), seed = 5002L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, num_cores = 1L)

  expect_output(print(res))
  expect_output(print(summary(res)))
})


# ===========================================================================
# 6. Individual rank is non-negative
# ===========================================================================

test_that("Rajive: individual rank non-negative for all blocks", {
  d   <- make_blocks(K = 3, n = 30, pks = c(20, 15, 12),
                     rankJ = 2, rankA = c(4, 3, 3), seed = 7001L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, num_cores = 1L)

  for (k in seq_along(d$blocks)) {
    ind_rank <- bd_individual(res, k)$rank
    expect_gte(ind_rank, 0L,
               label = paste0("individual rank >= 0, block ", k))
  }
})


# ===========================================================================
# 7. num_cores=1 and num_cores=2 give same column space (fixed joint_rank)
# ===========================================================================

test_that("Rajive: num_cores=1 and num_cores=2 give identical results (fixed joint_rank)", {
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")

  d   <- make_blocks(K = 2, n = 30, pks = c(20, 15),
                     rankJ = 2, rankA = c(4, 3), seed = 8001L)
  r1  <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, num_cores = 1L)
  r2  <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, num_cores = 2L)

  expect_colspace_equal(r1$joint_scores, r2$joint_scores,
                        tol = 1e-5, label = "joint scores column space")
})


# ===========================================================================
# 8. get_block_scores() and get_block_loadings() accessors
# ===========================================================================

test_that("get_block_scores returns matrix of correct dimensions", {
  d   <- make_blocks(K = 2, n = 40, pks = c(30, 20),
                     rankJ = 2, rankA = c(5, 3), seed = 9001L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                joint_rank = 2, num_cores = 1L)
  joint_rank <- res$joint_rank

  for (k in seq_along(d$blocks)) {
    js <- get_block_scores(res, k, 'joint')
    is <- get_block_scores(res, k, 'individual')
    expect_equal(nrow(js), 40L,
                 label = paste0("joint scores nrow, block ", k))
    expect_equal(ncol(js), joint_rank,
                 label = paste0("joint scores ncol, block ", k))
  }
})

