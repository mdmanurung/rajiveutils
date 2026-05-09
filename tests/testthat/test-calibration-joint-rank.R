# Slow calibration test for the Rajive() joint-rank selection rule.
#
# W-M8 (calibration half): Under the global null (independent Gaussian
# blocks, no joint structure), the Type-I rate for "joint_rank > 0" should
# be bounded, documenting the actual empirical performance of the
# max(wedin, rand_dir | perm) heuristic.
#
# Run with:
#   RAJIVE_RUN_SLOW=1 conda run -n R4_51 Rscript -e \
#     "devtools::test_file('tests/testthat/test-calibration-joint-rank.R')"
#
# Reference: Lock EF et al. (2013) Ann. Appl. Statist. 7(1):523-542.
# See also: R/Rajive.R @section "Joint-rank threshold (heuristic, no formal FWER/FDR)".

# ---------------------------------------------------------------------------
# W-M8: Type-I rate for joint_rank > 0 under H0
# ---------------------------------------------------------------------------

test_that("joint_rank Type-I rate under H0 is bounded at <= 0.10", {
  skip_if_not_slow()

  # B=200 replications, K=3 blocks.
  # Binomial sd at p=0.05, B=200: sqrt(0.05*0.95/200) ≈ 0.0154.
  # Upper 3-sigma bound ≈ 0.096; we assert <= 0.10 for a safety margin.
  with_lecuyer_seed(2026, {
    B <- 200L
    ranks <- replicate(B, {
      blk <- lapply(c(50L, 40L, 30L), function(p)
        matrix(stats::rnorm(60 * p), 60, p))
      Rajive(blk, c(6L, 5L, 4L),
             n_wedin_samples     = 200L,
             n_rand_dir_samples  = 200L,
             num_cores           = 1L)$joint_rank
    })

    type1_rate <- mean(ranks > 0)
    cat(sprintf("\n  [W-M8 calibration] joint_rank > 0 rate = %.4f (B=%d)\n",
                type1_rate, B))

    expect_lte(type1_rate, 0.10,
               label = paste0("Type-I rate = ", signif(type1_rate, 3),
                              "; must be <= 0.10"))
  })
})

# ---------------------------------------------------------------------------
# W-M8: Power — joint_rank >= true rank under strong signal
# ---------------------------------------------------------------------------

test_that("joint_rank recovers true rank >= 90% of the time under signal", {
  skip_if_not_slow()

  # B=100 reps; binomial sd at p=0.90, B=100 ≈ 0.030; lower 3-sigma ≈ 0.81.
  # We assert >= 0.80 for a safety margin.
  with_lecuyer_seed(4242, {
    B <- 100L
    correct <- replicate(B, {
      Y   <- ajive.data.sim(K = 3, rankJ = 2, rankA = c(5, 4, 3),
                            n = 60, pks = c(50, 40, 30), dist.type = 1)
      fit <- Rajive(Y$sim_data, c(5, 4, 3),
                   n_wedin_samples    = 200L,
                   n_rand_dir_samples = 200L,
                   num_cores          = 1L)
      fit$joint_rank >= 2L
    })

    power <- mean(correct)
    cat(sprintf("  [W-M8 power] joint_rank >= 2 rate = %.4f (B=%d)\n",
                power, B))

    expect_gte(power, 0.80,
               label = paste0("Power = ", signif(power, 3),
                              "; must be >= 0.80"))
  })
})
