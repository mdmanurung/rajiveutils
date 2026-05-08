library(testthat)
library(rajiveplus)

# ---------------------------------------------------------------------------
# Helpers shared across tests
# ---------------------------------------------------------------------------

make_simple_ajive <- function(n = 40, d1 = 30, d2 = 20, joint_rank = 2,
                              seed = 1L) {
  set.seed(seed)
  pks <- c(d1, d2)
  Y   <- ajive.data.sim(K = 2, rankJ = joint_rank, rankA = c(4, 3),
                        n = n, pks = pks, dist.type = 1)
  data_blocks          <- Y$sim_data
  initial_signal_ranks <- c(4, 3)
  ajive_out <- Rajive(data_blocks, initial_signal_ranks)
  list(ajive_out = ajive_out, blocks = data_blocks)
}

# ---------------------------------------------------------------------------
# 1. ols_f_stat_matrix correctness
# ---------------------------------------------------------------------------

test_that("ols_f_stat_matrix matches lm() F-statistic", {
  set.seed(7)
  n  <- 50
  x  <- rnorm(n)
  y1 <- 2 * x + rnorm(n)
  y2 <- rnorm(n)                      # unrelated feature

  Y_t <- rbind(y1, y2)               # 2 x n

  f_mat <- rajiveplus:::ols_f_stat_matrix(Y_t, x)

  # Reference from lm()
  f_lm1 <- summary(lm(y1 ~ x))$fstatistic[1]
  f_lm2 <- summary(lm(y2 ~ x))$fstatistic[1]

  expect_equal(as.numeric(f_mat[1]), as.numeric(f_lm1), tolerance = 1e-8)
  expect_equal(as.numeric(f_mat[2]), as.numeric(f_lm2), tolerance = 1e-8)
})

test_that("ols_f_stat_matrix returns NA for constant features", {
  n   <- 30
  x   <- rnorm(n)
  Y_t <- rbind(rep(5, n), rnorm(n))
  f   <- rajiveplus:::ols_f_stat_matrix(Y_t, x)
  expect_true(is.na(f[1]))
  expect_false(is.na(f[2]))
})

# ---------------------------------------------------------------------------
# 2. compute_empirical_pvalues
# ---------------------------------------------------------------------------

test_that("compute_empirical_pvalues returns 1 for NA observed F-stat", {
  f_obs  <- c(NA_real_, 5.0)
  f_null <- matrix(c(1, 2, 3, 4, 1, 2, 3, 4), nrow = 2, ncol = 4)
  p <- rajiveplus:::compute_empirical_pvalues(f_obs, f_null)
  expect_equal(p[1], 1)
})

test_that("compute_empirical_pvalues: p-value is proportion of nulls exceeding obs", {
  set.seed(1)
  # observed F >> null => p should be ~0
  f_obs  <- c(100)
  f_null <- matrix(runif(100, 0, 5), nrow = 1, ncol = 100)
  p <- rajiveplus:::compute_empirical_pvalues(f_obs, f_null)
  expect_lt(p[1], 0.05)
})

test_that("compute_empirical_pvalues: p-value is ~1 when obs F << nulls", {
  set.seed(2)
  f_obs  <- c(0.001)
  f_null <- matrix(runif(100, 10, 20), nrow = 1, ncol = 100)
  p <- rajiveplus:::compute_empirical_pvalues(f_obs, f_null)
  expect_equal(p[1], 1)
})

# ---------------------------------------------------------------------------
# 3. jackstraw_rajive structure
# ---------------------------------------------------------------------------

test_that("jackstraw_rajive returns correct structure", {
  res <- make_simple_ajive()
  js  <- jackstraw_rajive(res$ajive_out, res$blocks,
                          alpha = 0.05, n_null = 5, correction = "bonferroni")

  expect_s3_class(js, "jackstraw_rajive")

  joint_rank <- attr(js, "joint_rank")
  K          <- attr(js, "n_blocks")
  expect_equal(K, 2L)
  expect_equal(attr(js, "n_tests"), sum(vapply(res$blocks, ncol, integer(1L))) * joint_rank)

  # Block names
  expect_equal(names(js), c("block1", "block2"))

  for (k in seq_len(K)) {
    for (j in seq_len(joint_rank)) {
      comp <- js[[k]][[j]]
      # Expected fields
      expect_true(all(c("f_obs", "f_null", "p_values", "p_adj",
                        "significant", "significant_vars") %in% names(comp)))
      d_k <- ncol(res$blocks[[k]])
      expect_length(comp$f_obs,    d_k)
      expect_length(comp$p_values, d_k)
      expect_length(comp$p_adj,    d_k)
      expect_length(comp$significant, d_k)
      expect_equal(nrow(comp$f_null), d_k)
      expect_equal(ncol(comp$f_null), 5L)
    }
  }
})

# ---------------------------------------------------------------------------
# 4. Global correction behavior
# ---------------------------------------------------------------------------

test_that("Bonferroni correction is more conservative than no correction", {
  # Use the same seed before each call so the null distributions are identical.
  # Under the same empirical p-values, global Bonferroni is always
  # a stricter threshold than "none" (alpha), so significant_bonf is a
  # subset of significant_none.
  res <- make_simple_ajive(seed = 99L)

  set.seed(123L)
  js_bonf <- jackstraw_rajive(res$ajive_out, res$blocks,
                              alpha = 0.05, n_null = 10,
                              correction = "bonferroni")
  set.seed(123L)
  js_none <- jackstraw_rajive(res$ajive_out, res$blocks,
                              alpha = 0.05, n_null = 10,
                              correction = "none")

  for (k in 1:2) {
    for (j in seq_len(attr(js_bonf, "joint_rank"))) {
      # With identical p-values, Bonferroni can only be equal or more strict
      n_sig_bonf <- sum(js_bonf[[k]][[j]]$significant, na.rm = TRUE)
      n_sig_none <- sum(js_none[[k]][[j]]$significant, na.rm = TRUE)
      expect_lte(n_sig_bonf, n_sig_none)
    }
  }
})

test_that("BH adjustment is applied globally across blocks and components", {
  res <- make_simple_ajive(seed = 123L)

  set.seed(321L)
  js_bh <- jackstraw_rajive(res$ajive_out, res$blocks,
                            alpha = 0.05, n_null = 10,
                            correction = "BH")

  raw_p <- unlist(lapply(js_bh, function(block_result) {
    unlist(lapply(block_result, `[[`, "p_values"), use.names = FALSE)
  }), use.names = FALSE)
  expected_adj <- p.adjust(raw_p, method = "BH")
  observed_adj <- unlist(lapply(js_bh, function(block_result) {
    unlist(lapply(block_result, `[[`, "p_adj"), use.names = FALSE)
  }), use.names = FALSE)

  expect_equal(observed_adj, expected_adj)
})

test_that("default correction is BH", {
  res <- make_simple_ajive(seed = 11L)
  js <- jackstraw_rajive(res$ajive_out, res$blocks, n_null = 5)

  expect_identical(attr(js, "correction"), "BH")
})

# ---------------------------------------------------------------------------
# 5. Degenerate-feature guard (W-M10)
# ---------------------------------------------------------------------------

test_that("constant feature triggers degenerate-block error", {
  set.seed(42)
  n  <- 60
  pks <- c(20, 15)
  # Use simulated data with a true joint structure so the identifiability
  # filter does not drop every joint component (joint_rank must be > 0).
  Y <- ajive.data.sim(K = 2, rankJ = 1, rankA = c(3, 3), n = n,
                      pks = pks, dist.type = 1)
  blocks <- Y$sim_data
  # Inject a constant column into block 1
  blocks[[1]][, 1] <- 7.0          # constant feature

  expect_error(
    Rajive(blocks, initial_signal_ranks = c(3, 3)),
    class = "rajiveplus_degenerate_block"
  )
})

# ---------------------------------------------------------------------------
# 6. p-values for strongly correlated feature should be small
# ---------------------------------------------------------------------------

test_that("feature strongly correlated with joint scores gets a small p-value", {
  set.seed(5)
  n  <- 60
  pks <- c(50, 40)
  Y   <- ajive.data.sim(K = 2, rankJ = 1, rankA = c(3, 3), n = n,
                        pks = pks, dist.type = 1)
  blocks <- Y$sim_data
  ajive_out <- Rajive(blocks, initial_signal_ranks = c(3, 3))

  # Replace first feature of block 1 with a near-perfect copy of joint score
  joint_sc <- ajive_out$joint_scores[, 1]
  blocks[[1]][, 1] <- joint_sc + rnorm(n, sd = 0.001)

  js <- jackstraw_rajive(ajive_out, blocks,
                         alpha = 0.05, n_null = 20, correction = "none")

  # First feature should have a very small p-value (or at least rank 1st)
  p1 <- js[["block1"]][["comp1"]]$p_values
  expect_equal(which.min(p1), 1L)
})

# ---------------------------------------------------------------------------
# 7. get_significant_vars accessor
# ---------------------------------------------------------------------------

test_that("get_significant_vars returns the correct vector", {
  res <- make_simple_ajive(seed = 3L)
  js  <- jackstraw_rajive(res$ajive_out, res$blocks,
                          alpha = 0.05, n_null = 5, correction = "bonferroni")

  sv <- get_significant_vars(js, block = 1, component = 1)

  # Must equal the stored significant_vars in the nested list
  expect_identical(sv, js[["block1"]][["comp1"]]$significant_vars)
})

test_that("get_significant_vars raises error on bad block index", {
  res <- make_simple_ajive()
  js  <- jackstraw_rajive(res$ajive_out, res$blocks, n_null = 5)
  expect_error(get_significant_vars(js, block = 99, component = 1),
               regexp = "block")
})

# ---------------------------------------------------------------------------
# 8. Variable names are propagated when colnames are present
# ---------------------------------------------------------------------------

test_that("variable names are attached to results when colnames are present", {
  set.seed(10)
  n    <- 60
  # Use data with a true joint signal so joint_rank > 0 after the
  # (L2-norm) identifiability filter.
  Y <- ajive.data.sim(K = 2, rankJ = 1, rankA = c(2, 2), n = n,
                      pks = c(10, 8), dist.type = 1)
  blocks <- Y$sim_data
  colnames(blocks[[1]]) <- paste0("gene", seq_len(ncol(blocks[[1]])))
  colnames(blocks[[2]]) <- paste0("protein", seq_len(ncol(blocks[[2]])))

  ajive_out <- Rajive(blocks, initial_signal_ranks = c(2, 2))
  js <- jackstraw_rajive(ajive_out, blocks, n_null = 5, correction = "none")

  comp1 <- js[["block1"]][["comp1"]]
  expect_equal(names(comp1$f_obs),    colnames(blocks[[1]]))
  expect_equal(names(comp1$p_values), colnames(blocks[[1]]))
})

# ---------------------------------------------------------------------------
# 9. Input validation
# ---------------------------------------------------------------------------

test_that("jackstraw_rajive raises informative errors on bad inputs", {
  res <- make_simple_ajive()

  expect_error(jackstraw_rajive(res$ajive_out, list()),
               regexp = "non-empty")
  expect_error(jackstraw_rajive(res$ajive_out, res$blocks, alpha = 1.5),
               regexp = "alpha")
  expect_error(jackstraw_rajive(res$ajive_out, res$blocks, n_null = 0),
               regexp = "n_null")
})

# ---------------------------------------------------------------------------
# 10. print and summary do not error
# ---------------------------------------------------------------------------

test_that("print.jackstraw_rajive and summary.jackstraw_rajive run without error", {
  res <- make_simple_ajive()
  js  <- jackstraw_rajive(res$ajive_out, res$blocks, n_null = 5)
  expect_output(print(js))
  df <- summary(js)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("block", "component", "n_features", "n_significant",
                    "alpha", "correction") %in% names(df)))
})

# ---------------------------------------------------------------------------
# 11. Phipson--Smyth p-value correctness (regression test for audit fix #1)
# ---------------------------------------------------------------------------

test_that("compute_empirical_pvalues never returns 0 (Phipson--Smyth)", {
  set.seed(11)
  # Observed F-stats far above the null pool: would have given p = 0
  # under the old formula; with Phipson--Smyth they must be > 0.
  f_obs  <- c(1e6, 1e6, 1e6)
  f_null <- matrix(runif(3 * 50, 0, 1), nrow = 3, ncol = 50)
  p <- rajiveplus:::compute_empirical_pvalues(f_obs, f_null)
  expect_true(all(p > 0))
  # Lower bound = 1 / (1 + N_pool) = 1 / 151
  expect_equal(min(p), 1 / (1 + 3 * 50))
})

test_that("compute_empirical_pvalues pools nulls across features", {
  # With a single observed F-stat and all nulls in one feature row,
  # pooling across rows must use the full pool, not just that row.
  f_obs  <- 1.5
  # Place pool values across multiple rows; obs is below the pool max.
  f_null <- matrix(c(rep(NA_real_, 3), c(1.0, 1.2, 1.4, 1.6, 1.8, 2.0)),
                   nrow = 3, ncol = 3, byrow = FALSE)
  # Pool (after dropping NAs) = c(1.0, 1.2, 1.4, 1.6, 1.8, 2.0); N=6.
  # n_ge(1.5) = #{>= 1.5} = 3 (1.6, 1.8, 2.0). p = (1+3)/(1+6) = 4/7.
  p <- rajiveplus:::compute_empirical_pvalues(f_obs, f_null)
  expect_equal(p, 4 / 7)
})

test_that("jackstraw p-values are bounded away from 0 and well-spread", {
  # Regression test for audit fix #1: the old formula `mean(f_obs < nulls)`
  # could produce p_value = 0, which inflates FDR after BH adjustment.  With
  # Phipson--Smyth pooling, p-values are bounded below by 1 / (1 + N_pool)
  # and the empirical distribution is finer-grained than the per-feature
  # null pool would allow.
  set.seed(123)
  n  <- 60
  X1 <- matrix(rnorm(n * 80), n, 80)
  X2 <- matrix(rnorm(n * 60), n, 60)
  blocks <- list(X1, X2)
  ajive_out <- Rajive(blocks, initial_signal_ranks = c(3, 3), joint_rank = 2)
  js <- jackstraw_rajive(ajive_out, blocks,
                         alpha = 0.05, n_null = 20, correction = "none")

  # Expected lower bound per block-component: 1 / (1 + d_k * n_null).
  # For block 1: 1 / (1 + 80*20) = 1/1601.
  # For block 2: 1 / (1 + 60*20) = 1/1201.
  for (k in seq_along(js)) {
    d_k <- length(blocks[[k]][1, ])
    lb  <- 1 / (1 + d_k * 20)
    for (j in seq_len(attr(js, "joint_rank"))) {
      p <- js[[k]][[j]]$p_values
      p <- p[!is.na(p)]
      # No exact zeros (regression test for old bug).
      expect_true(all(p > 0))
      # Phipson--Smyth lower bound respected.
      expect_gte(min(p), lb - 1e-12)
      # P-values are well-spread (more than n_null+1 distinct values),
      # confirming pooling rather than per-feature CDF.
      expect_gt(length(unique(p)), 20L + 1L)
    }
  }
})

# ---------------------------------------------------------------------------
# 12. W-M1: compute_empirical_pvalues uses >= (exact-tie correctness)
# ---------------------------------------------------------------------------

test_that("compute_empirical_pvalues uses >= not > (W-M1 regression)", {
  # With pool = {1,1,2,2,3,3} and f_obs = {1, 2, 3}:
  #   #{pool >= 1} = 6, #{pool >= 2} = 4, #{pool >= 3} = 2.
  # Under the old findInterval(left.open=FALSE): #{pool > x} = {4, 2, 0},
  # so p would be {5/7, 3/7, 1/7} -- wrong for tied values.
  # After fix (left.open=TRUE): #{pool >= x} is correct.
  f_obs  <- c(1, 2, 3)
  f_null <- matrix(c(1, 1, 2, 2, 3, 3), nrow = 3)
  # Pool = c(1,1,2,2,3,3), N = 6.
  expected <- (1 + c(6, 4, 2)) / (1 + 6)
  got <- rajiveplus:::compute_empirical_pvalues(f_obs, f_null)
  expect_equal(got, expected)
})

test_that("compute_empirical_pvalues >= formula is identity for continuous pool", {
  # For a continuous distribution, #{pool >= x} ≈ #{pool > x} for any given x,
  # but we verify the formula gives consistent results with reference.
  set.seed(99)
  f_obs  <- c(1.5, 3.0, 5.0)
  pool   <- sort(runif(999, 0, 10))
  f_null <- matrix(pool, nrow = length(f_obs), ncol = length(pool) / length(f_obs),
                   byrow = TRUE)
  # Reference: brute-force #{pool_full >= f_obs[i]} / (1 + length(pool_full))
  pool_full <- as.numeric(f_null)
  expected <- (1 + vapply(f_obs, function(x) sum(pool_full >= x), integer(1L))) /
              (1 + length(pool_full))
  got <- rajiveplus:::compute_empirical_pvalues(f_obs, f_null)
  expect_equal(got, expected, tolerance = 0)
})
