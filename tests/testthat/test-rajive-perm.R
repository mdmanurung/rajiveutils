library(testthat)
library(rajiveplus)

# ---------------------------------------------------------------------------
# Tests for n_perm_samples (permutation-based joint rank threshold)
#
# Strategy: small simulated blocks for speed; modest n_perm_samples (50)
# to keep runtime low while still exercising the full code path.
# ---------------------------------------------------------------------------

make_blocks_perm <- function(K, n, pks, rankJ, rankA, seed) {
  set.seed(seed)
  Y <- ajive.data.sim(K = K, rankJ = rankJ, rankA = rankA,
                      n = n, pks = pks, dist.type = 1)
  list(blocks = Y$sim_data, initial_signal_ranks = rankA)
}


# ===========================================================================
# 1. Backwards compatibility: NA default produces unchanged behaviour
# ===========================================================================

test_that("n_perm_samples=NA (default) preserves original rand_dir behaviour", {
  d <- make_blocks_perm(K = 2, n = 30, pks = c(20, 15),
                        rankJ = 2, rankA = c(4, 3), seed = 7001L)

  set.seed(42L)
  res_default <- Rajive(d$blocks, d$initial_signal_ranks,
                        n_wedin_samples = 30L, n_rand_dir_samples = 30L,
                        num_cores = 1L)

  set.seed(42L)
  res_explicit <- Rajive(d$blocks, d$initial_signal_ranks,
                         n_wedin_samples = 30L, n_rand_dir_samples = 30L,
                         n_perm_samples = NA,
                         num_cores = 1L)

  # Same rand_dir path -> same diagnostic structure
  expect_true(!is.null(res_default$joint_rank_sel$rand_dir))
  expect_true(!is.null(res_explicit$joint_rank_sel$rand_dir))
  expect_null(res_default$joint_rank_sel$perm)
  expect_null(res_explicit$joint_rank_sel$perm)
})


# ===========================================================================
# 2. n_perm_samples non-NA populates rank_sel_results$perm
# ===========================================================================

test_that("n_perm_samples populates joint_rank_sel$perm with expected fields", {
  d <- make_blocks_perm(K = 2, n = 30, pks = c(20, 15),
                        rankJ = 2, rankA = c(4, 3), seed = 7002L)

  set.seed(11L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                n_wedin_samples = 30L,
                n_perm_samples = 50L,
                num_cores = 1L)

  expect_true(!is.null(res$joint_rank_sel$perm))
  expect_named(res$joint_rank_sel$perm,
               c("perm_samples", "perm_svsq_threshold"))
  expect_length(res$joint_rank_sel$perm$perm_samples, 50L)
  expect_true(is.numeric(res$joint_rank_sel$perm$perm_svsq_threshold))
  expect_true(res$joint_rank_sel$perm$perm_svsq_threshold > 0)
})


# ===========================================================================
# 3. perm replaces rand_dir (mutually exclusive)
# ===========================================================================

test_that("when n_perm_samples is set, rand_dir is NOT computed", {
  d <- make_blocks_perm(K = 2, n = 30, pks = c(20, 15),
                        rankJ = 2, rankA = c(4, 3), seed = 7003L)

  set.seed(11L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                n_wedin_samples = 30L,
                n_rand_dir_samples = 30L,
                n_perm_samples = 50L,
                num_cores = 1L)

  expect_true(!is.null(res$joint_rank_sel$perm))
  expect_null(res$joint_rank_sel$rand_dir)
})


# ===========================================================================
# 4. Threshold combination uses max(wedin, perm)
# ===========================================================================

test_that("overall_sv_sq_threshold equals max(wedin, perm)", {
  d <- make_blocks_perm(K = 2, n = 30, pks = c(20, 15),
                        rankJ = 2, rankA = c(4, 3), seed = 7004L)

  set.seed(11L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                n_wedin_samples = 30L,
                n_perm_samples = 50L,
                num_cores = 1L)

  jrs <- res$joint_rank_sel
  expected <- max(jrs$wedin$wedin_svsq_threshold,
                  jrs$perm$perm_svsq_threshold,
                  na.rm = TRUE)
  expect_equal(as.numeric(jrs$overall_sv_sq_threshold),
               as.numeric(expected))
})


# ===========================================================================
# 5. Recovers known joint rank using perm-only threshold (no rand_dir, no wedin)
# ===========================================================================

test_that("perm-only path recovers known joint rank on strong-signal sim", {
  # Use a stronger-signal simulation so the permutation null is clearly
  # separated from the observed joint singular values.
  d <- make_blocks_perm(K = 2, n = 80, pks = c(60, 50),
                        rankJ = 2, rankA = c(4, 3), seed = 7005L)

  set.seed(11L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                n_wedin_samples = NA,
                n_rand_dir_samples = NA,
                n_perm_samples = 100L,
                num_cores = 1L)

  expect_true(!is.null(res$joint_rank_sel$perm))
  expect_null(res$joint_rank_sel$wedin)
  expect_null(res$joint_rank_sel$rand_dir)
  # Permutation is a conservative non-parametric null; verify it recovers
  # at least the strongest joint component and does not exceed the truth.
  expect_gte(res$joint_rank, 1L)
  expect_lte(res$joint_rank, 2L)
})


# ===========================================================================
# 6. Error guard: all four NA -> error
# ===========================================================================

test_that("error when all rank-selection arguments are NA", {
  d <- make_blocks_perm(K = 2, n = 30, pks = c(20, 15),
                        rankJ = 2, rankA = c(4, 3), seed = 7006L)

  expect_error(
    Rajive(d$blocks, d$initial_signal_ranks,
           n_wedin_samples = NA,
           n_rand_dir_samples = NA,
           n_perm_samples = NA,
           joint_rank = NA,
           num_cores = 1L),
    regexp = "n_perm_samples"
  )
})


# ===========================================================================
# 7. Reproducibility: identical seed -> identical perm samples
# ===========================================================================

test_that("perm samples are reproducible under set.seed", {
  d <- make_blocks_perm(K = 2, n = 30, pks = c(20, 15),
                        rankJ = 2, rankA = c(4, 3), seed = 7007L)

  set.seed(123L)
  r1 <- Rajive(d$blocks, d$initial_signal_ranks,
               n_wedin_samples = NA,
               n_perm_samples = 30L,
               num_cores = 1L)

  set.seed(123L)
  r2 <- Rajive(d$blocks, d$initial_signal_ranks,
               n_wedin_samples = NA,
               n_perm_samples = 30L,
               num_cores = 1L)

  expect_equal(r1$joint_rank_sel$perm$perm_samples,
               r2$joint_rank_sel$perm$perm_samples)
})


# ===========================================================================
# 8. Visualization payload exposes has_perm and perm_cutoff
# ===========================================================================

test_that("extract_components(what='rank_diagnostics') surfaces perm fields", {
  d <- make_blocks_perm(K = 2, n = 30, pks = c(20, 15),
                        rankJ = 2, rankA = c(4, 3), seed = 7008L)

  set.seed(7L)
  res <- Rajive(d$blocks, d$initial_signal_ranks,
                n_wedin_samples = 30L,
                n_perm_samples = 30L,
                num_cores = 1L)

  diag <- extract_components(res, what = "rank_diagnostics",
                             format = "wide")
  expect_true(diag$has_perm)
  expect_false(diag$has_random)
  expect_true(is.numeric(diag$perm_cutoff))
  expect_length(diag$perm_samples, 30L)
  expect_match(diag$cutoff_rule, "perm", fixed = TRUE)
})
