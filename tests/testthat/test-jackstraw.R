library(testthat)
library(RaJIVE)

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

  f_mat <- RaJIVE:::ols_f_stat_matrix(Y_t, x)

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
  f   <- RaJIVE:::ols_f_stat_matrix(Y_t, x)
  expect_true(is.na(f[1]))
  expect_false(is.na(f[2]))
})

# ---------------------------------------------------------------------------
# 2. compute_empirical_pvalues
# ---------------------------------------------------------------------------

test_that("compute_empirical_pvalues returns 1 for NA observed F-stat", {
  f_obs  <- c(NA_real_, 5.0)
  f_null <- matrix(c(1, 2, 3, 4, 1, 2, 3, 4), nrow = 2, ncol = 4)
  p <- RaJIVE:::compute_empirical_pvalues(f_obs, f_null)
  expect_equal(p[1], 1)
})

test_that("compute_empirical_pvalues: p-value is proportion of nulls exceeding obs", {
  set.seed(1)
  # observed F >> null => p should be ~0
  f_obs  <- c(100)
  f_null <- matrix(runif(100, 0, 5), nrow = 1, ncol = 100)
  p <- RaJIVE:::compute_empirical_pvalues(f_obs, f_null)
  expect_lt(p[1], 0.05)
})

test_that("compute_empirical_pvalues: p-value is ~1 when obs F << nulls", {
  set.seed(2)
  f_obs  <- c(0.001)
  f_null <- matrix(runif(100, 10, 20), nrow = 1, ncol = 100)
  p <- RaJIVE:::compute_empirical_pvalues(f_obs, f_null)
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
# 4. Bonferroni gives fewer or equal significant variables than "none"
# ---------------------------------------------------------------------------

test_that("Bonferroni correction is more conservative than no correction", {
  # Use the same seed before each call so the null distributions are identical.
  # Under the same empirical p-values, Bonferroni (alpha/(d*jr)) is always
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

# ---------------------------------------------------------------------------
# 5. Constant feature receives p-value = 1 and is not significant
# ---------------------------------------------------------------------------

test_that("constant feature is not flagged as significant", {
  set.seed(42)
  n  <- 40
  d  <- 20
  jr <- 1L

  # Inject a constant column into block 1
  X1    <- matrix(rnorm(n * d), n, d)
  X1[, 1] <- 7.0          # constant feature

  X2    <- matrix(rnorm(n * 15), n, 15)

  # Build a minimal ajive_output by hand using RaJIVE internals
  # (run actual Rajive so we get a valid joint_scores)
  blocks <- list(X1, X2)
  ajive_out <- Rajive(blocks, initial_signal_ranks = c(3, 3))

  js <- jackstraw_rajive(ajive_out, blocks,
                         alpha = 0.05, n_null = 5, correction = "bonferroni")

  # Feature 1 in block 1 is constant -> p-value == 1 -> not significant
  for (j in seq_len(attr(js, "joint_rank"))) {
    expect_equal(js[["block1"]][[j]]$p_values[1], 1)
    expect_false(isTRUE(js[["block1"]][[j]]$significant[1]))
  }
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
  n    <- 40
  X1   <- matrix(rnorm(n * 10), n, 10)
  X2   <- matrix(rnorm(n * 8),  n, 8)
  colnames(X1) <- paste0("gene", 1:10)
  colnames(X2) <- paste0("protein", 1:8)

  blocks <- list(X1, X2)
  ajive_out <- Rajive(blocks, initial_signal_ranks = c(2, 2))
  js <- jackstraw_rajive(ajive_out, blocks, n_null = 5, correction = "none")

  comp1 <- js[["block1"]][["comp1"]]
  expect_equal(names(comp1$f_obs),    colnames(X1))
  expect_equal(names(comp1$p_values), colnames(X1))
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
