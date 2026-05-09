# Slow calibration tests for jackstraw_rajive().
#
# These tests verify that jackstraw p-values are Uniform(0,1) under the
# global null hypothesis (independent Gaussian blocks with no joint structure).
#
# Run with:
#   RAJIVE_RUN_SLOW=1 conda run -n R4_51 Rscript -e \
#     "devtools::test_file('tests/testthat/test-calibration-jackstraw.R')"
#
# Reference: Yang X et al. (2021) arXiv:2109.12272 (YHHM 2021).
# Reference: Phipson B, Smyth GK (2010) Stat. Appl. Genet. Mol. Biol. 9(1).

# ---------------------------------------------------------------------------
# W-M4: Null uniformity
# ---------------------------------------------------------------------------

test_that("jackstraw_rajive p-values are uniform under H0", {
  skip_if_not_slow()

  # Parameters: B=200 outer reps, n_null=50 per feature.
  # n_features = 40+30 = 70 per outer rep -> N = 70*200 = 14000 p-values.
  # Binomial sd at p=0.05: sqrt(0.05*0.95/14000) ≈ 0.0018.
  # 3-sigma CI ≈ [0.044, 0.056]; we use [0.040, 0.060] as safety margin.
  with_lecuyer_seed(2026, {
    B  <- 200L
    ps <- unlist(replicate(B, {
      blk <- list(
        matrix(stats::rnorm(60 * 40), 60, 40),
        matrix(stats::rnorm(60 * 30), 60, 30)
      )
      fit <- Rajive(blk, c(5, 4), joint_rank = 1L, num_cores = 1L)
      js  <- jackstraw_rajive(fit, blk, n_null = 50L, correction = "none")
      unlist(lapply(js, function(b) lapply(b, `[[`, "p_values")),
             use.names = FALSE)
    }, simplify = FALSE))

    ks_p <- stats::ks.test(ps, "punif")$p.value
    alpha_rate <- mean(ps <= 0.05)

    expect_gt(ks_p, 0.01,
              label = paste0("KS p-value = ", signif(ks_p, 3),
                             "; p-values should look uniform under H0"))
    expect_gte(alpha_rate, 0.040,
               label = paste0("Type-I rate = ", signif(alpha_rate, 3)))
    expect_lte(alpha_rate, 0.060,
               label = paste0("Type-I rate = ", signif(alpha_rate, 3)))
  })
})

# ---------------------------------------------------------------------------
# W-M4: Power property — significant features under strong signal
# ---------------------------------------------------------------------------

test_that("jackstraw_rajive detects >= 80% of true-loading features under signal", {
  skip_if_not_slow()

  # B=50 reps; binomial sd at p=0.9, B=50 ≈ 0.042; lower 3-sigma ≈ 0.77.
  # We assert >= 0.80 for a safety margin.
  with_lecuyer_seed(9999, {
    B         <- 50L
    detect    <- replicate(B, {
      Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4),
                            n = 60, pks = c(80, 60), dist.type = 1)
      fit <- Rajive(Y$sim_data, c(5, 4), num_cores = 1L)
      js  <- jackstraw_rajive(fit, Y$sim_data,
                              n_null = 50L, correction = "none")
      # fraction of block1 features called significant at alpha=0.05
      mean(js$block1$comp1$p_values <= 0.05, na.rm = TRUE)
    })

    mean_detect <- mean(detect)
    expect_gte(mean_detect, 0.80,
               label = paste0("Mean detection rate = ", signif(mean_detect, 3)))
  })
})
