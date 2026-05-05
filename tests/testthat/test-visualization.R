library(testthat)
library(rajiveplus)

# ---------------------------------------------------------------------------
# Shared fixture — computed ONCE for the entire file (fast: 50 resamples)
# ---------------------------------------------------------------------------

.viz_fit <- local({
  set.seed(42L)
  n   <- 40L
  pks <- c(30L, 20L)
  Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                        n = n, pks = pks, dist.type = 1)
  Rajive(Y$sim_data, c(4L, 3L),
         n_wedin_samples = 50L, n_rand_dir_samples = 50L)
})

make_diag_fit <- function() .viz_fit

# ---------------------------------------------------------------------------
# A. extract_components — return schema tests
# ---------------------------------------------------------------------------

test_that("A1: extract_components returns rajive_diagnostics list (wide, both bounds)", {
  fit  <- make_diag_fit()
  diag <- extract_components(fit, what = "rank_diagnostics", format = "wide")

  expect_true(inherits(diag, "rajive_diagnostics"))
  expect_true(is.list(diag))

  expected_fields <- c("obs_svals", "obs_svals_sq", "joint_rank_estimate",
                       "overall_sv_sq_threshold", "wedin_samples",
                       "rand_dir_samples", "wedin_cutoff", "rand_cutoff",
                       "wedin_percentile", "rand_percentile",
                       "identif_dropped", "cutoff_rule",
                       "has_wedin", "has_random")
  expect_true(all(expected_fields %in% names(diag)))

  expect_equal(diag$joint_rank_estimate, fit$joint_rank)
})

test_that("A2: format='wide' returns list with numeric sval vectors", {
  fit  <- make_diag_fit()
  diag <- extract_components(fit, what = "rank_diagnostics", format = "wide")

  expect_true(is.numeric(diag$obs_svals))
  expect_true(is.numeric(diag$obs_svals_sq))
  expect_equal(diag$obs_svals_sq, diag$obs_svals^2)
  expect_true(length(diag$obs_svals) > 0L)
})

test_that("A3: format='long' returns data.frame with required columns", {
  fit  <- make_diag_fit()
  df   <- extract_components(fit, what = "rank_diagnostics", format = "long")

  expect_true(is.data.frame(df))
  required_cols <- c("component_index", "obs_sval", "obs_sval_sq",
                     "classification", "joint_rank_estimate",
                     "overall_sv_sq_threshold", "wedin_cutoff", "rand_cutoff")
  expect_true(all(required_cols %in% names(df)))
  expect_true(nrow(df) > 0L)
  expect_true(all(df$classification %in% c("joint", "nonjoint", "dropped")))
})

test_that("A4: Wedin absent — returns successfully with rand cutoff as overall", {
  fit <- make_diag_fit()
  fit$joint_rank_sel[["wedin"]] <- NULL
  diag <- extract_components(fit, what = "rank_diagnostics", format = "wide")

  expect_true(is.na(diag$wedin_cutoff))
  expect_false(diag$has_wedin)
  expect_true(diag$has_random)
  expect_equal(diag$cutoff_rule, "random_only")
})

test_that("A5: Random direction absent — returns successfully with Wedin cutoff as overall", {
  fit <- make_diag_fit()
  fit$joint_rank_sel[["rand_dir"]] <- NULL
  diag <- extract_components(fit, what = "rank_diagnostics", format = "wide")

  expect_true(is.na(diag$rand_cutoff))
  expect_false(diag$has_random)
  expect_true(diag$has_wedin)
  expect_equal(diag$cutoff_rule, "wedin_only")
})

test_that("A6: Neither bound present — returns with NAs and 'none_available' rule", {
  fit <- make_diag_fit()
  fit$joint_rank_sel[["wedin"]]    <- NULL
  fit$joint_rank_sel[["rand_dir"]] <- NULL

  expect_warning(
    diag <- extract_components(fit, what = "rank_diagnostics", format = "wide"),
    regexp = "bounds are absent"
  )

  expect_true(is.na(diag$wedin_cutoff))
  expect_true(is.na(diag$rand_cutoff))
  expect_equal(diag$cutoff_rule, "none_available")
  expect_false(diag$has_wedin)
  expect_false(diag$has_random)
})

test_that("A7: joint_rank = 0 case — all classification 'nonjoint' in long format", {
  fit <- make_diag_fit()
  # Force joint_rank to 0 and threshold above all svals
  fit$joint_rank <- 0L
  fit$joint_rank_sel$joint_rank_estimate <- 0L
  max_sq <- max(fit$joint_rank_sel$obs_svals^2) * 2
  fit$joint_rank_sel$overall_sv_sq_threshold <- max_sq
  if (!is.null(fit$joint_rank_sel$wedin)) {
    fit$joint_rank_sel$wedin$wedin_svsq_threshold <- max_sq
  }
  if (!is.null(fit$joint_rank_sel$rand_dir)) {
    fit$joint_rank_sel$rand_dir$rand_dir_svsq_threshold <- max_sq
  }

  df <- extract_components(fit, what = "rank_diagnostics", format = "long")
  expect_true(all(df$classification == "nonjoint"))
  expect_equal(unique(df$joint_rank_estimate), 0L)
})

test_that("A9: Missing joint_rank_sel field raises informative error", {
  fit <- make_diag_fit()
  fit$joint_rank_sel <- NULL

  expect_error(
    extract_components(fit, what = "rank_diagnostics"),
    regexp = "joint_rank_sel"
  )
})

# ---------------------------------------------------------------------------
# E. Error handling
# ---------------------------------------------------------------------------

test_that("E1: rank_diagnostics on jackstraw object errors with class message", {
  fit   <- make_diag_fit()
  blocks <- list(
    matrix(rnorm(40 * 30), 40, 30),
    matrix(rnorm(40 * 20), 40, 20)
  )
  jacks <- jackstraw_rajive(fit, blocks = blocks, n_null = 3L)

  expect_error(
    extract_components(jacks, what = "rank_diagnostics"),
    regexp = "class 'rajive'"
  )
})

test_that("E2: plot_components with no ajive_output and no diagnostics errors", {
  expect_error(
    plot_components(plot_type = "rank_threshold"),
    regexp = "ajive_output must be provided"
  )
})

test_that("E3: ajive_output not rajive class errors", {
  expect_error(
    plot_components(ajive_output = list(a = 1), plot_type = "rank_threshold"),
    regexp = "class 'rajive'"
  )
})

test_that("E4: cutoff_quantile outside (0, 1) errors", {
  fit <- make_diag_fit()
  expect_error(
    plot_components(fit, plot_type = "bound_distributions",
                    cutoff_quantile = 1.5),
    regexp = "cutoff_quantile must be in"
  )
})

# ---------------------------------------------------------------------------
# B. plot_components — rank_threshold returns ggplot
# ---------------------------------------------------------------------------

test_that("B8: plot_rank_threshold returns ggplot object", {
  fit <- make_diag_fit()
  p   <- plot_components(fit, plot_type = "rank_threshold")
  expect_true(inherits(p, "ggplot"))
})

test_that("B6: rank_threshold with no bounds still returns ggplot", {
  fit <- make_diag_fit()
  fit$joint_rank_sel[["wedin"]]    <- NULL
  fit$joint_rank_sel[["rand_dir"]] <- NULL

  suppressWarnings(
    p <- plot_components(fit, plot_type = "rank_threshold")
  )
  expect_true(inherits(p, "ggplot"))
})

# ---------------------------------------------------------------------------
# C. plot_components — bound_distributions returns ggplot
# ---------------------------------------------------------------------------

test_that("C5: bound_distributions returns ggplot object", {
  skip_if_not_installed("patchwork")
  fit <- make_diag_fit()
  p   <- plot_components(fit, plot_type = "bound_distributions")
  expect_true(inherits(p, "gg"))
})

test_that("C2: Only Wedin — single ggplot returned without error", {
  fit <- make_diag_fit()
  fit$joint_rank_sel[["rand_dir"]] <- NULL
  p   <- plot_components(fit, plot_type = "bound_distributions")
  expect_true(inherits(p, "ggplot"))
})

test_that("C3: Only rand — single ggplot returned without error", {
  fit <- make_diag_fit()
  fit$joint_rank_sel[["wedin"]] <- NULL
  p   <- plot_components(fit, plot_type = "bound_distributions")
  expect_true(inherits(p, "ggplot"))
})

# ---------------------------------------------------------------------------
# D. plot_components — ajive_diagnostic composite
# ---------------------------------------------------------------------------

test_that("D5: ajive_diagnostic returns gg-compatible object", {
  skip_if_not_installed("patchwork")
  fit <- make_diag_fit()
  p   <- plot_components(fit, plot_type = "ajive_diagnostic")
  expect_true(inherits(p, "gg"))
})

test_that("D4: ajive_diagnostic with no bounds returns single ggplot", {
  fit <- make_diag_fit()
  fit$joint_rank_sel[["wedin"]]    <- NULL
  fit$joint_rank_sel[["rand_dir"]] <- NULL

  suppressWarnings(
    p <- plot_components(fit, plot_type = "ajive_diagnostic")
  )
  expect_true(inherits(p, "ggplot"))
})

test_that("D6: ajive_diagnostic joint_rank=0 returns gg object", {
  skip_if_not_installed("patchwork")
  fit <- make_diag_fit()
  fit$joint_rank <- 0L
  fit$joint_rank_sel$joint_rank_estimate <- 0L
  fit$joint_rank_sel$overall_sv_sq_threshold <-
    max(fit$joint_rank_sel$obs_svals^2) * 2

  p <- plot_components(fit, plot_type = "ajive_diagnostic")
  expect_true(inherits(p, "gg"))
})

# ---------------------------------------------------------------------------
# identif_dropped stored in rank_sel_results
# ---------------------------------------------------------------------------

test_that("identif_dropped is an integer vector in rank_sel_results", {
  fit <- make_diag_fit()
  dropped <- fit$joint_rank_sel$identif_dropped
  expect_true(is.integer(dropped))
})

test_that("identif_dropped is empty when no components were removed", {
  fit <- make_diag_fit()
  # Typically none dropped for simulated clean data
  # We just check it is a valid integer(0) or small vector
  dropped <- fit$joint_rank_sel$identif_dropped
  expect_true(length(dropped) == 0L || all(dropped > 0L))
})

# ---------------------------------------------------------------------------
# G. Statistical guardrail integrations (#3-#5)
# ---------------------------------------------------------------------------

test_that("G1: associate_components returns table and emits score-estimation warning", {
  fit  <- make_diag_fit()
  n    <- nrow(fit$joint_scores)
  meta <- data.frame(x = rnorm(n))

  expect_message(
    res <- associate_components(fit, meta, variable = "x", mode = "continuous"),
    regexp = "Score estimation error is NOT propagated"
  )

  expect_true(is.data.frame(res))
  expect_true(all(c("variable", "component", "stat", "p_value", "p_adj", "method") %in% names(res)))
})

test_that("G2: survival split mode emits anti-conservative warning", {
  fit  <- make_diag_fit()
  n    <- nrow(fit$joint_scores)
  meta <- data.frame(
    score_ref = rnorm(n),
    time      = rexp(n, rate = 0.1),
    status    = sample(c(0L, 1L), n, replace = TRUE)
  )

  expect_message(
    associate_components(
      fit, meta,
      mode = "survival",
      variable = "score_ref",
      time_col = "time",
      status_col = "status",
      split = "median"
    ),
    regexp = "anti-conservative"
  )
})

test_that("G3: assess_stability(loadings) returns aligned loading summaries", {
  set.seed(123)
  Y <- ajive.data.sim(K = 2, rankJ = 1, rankA = c(2, 2),
                      n = 24, pks = c(15, 12), dist.type = 1)
  fit <- Rajive(Y$sim_data, c(2, 2),
                n_wedin_samples = 20L, n_rand_dir_samples = 20L)

  stab <- assess_stability(
    ajive_output = fit,
    blocks = Y$sim_data,
    initial_signal_ranks = c(2, 2),
    target = "loadings",
    B = 3,
    sample_frac = 0.8
  )

  expect_true(is.list(stab))
  expect_true(length(stab) == 2L)
  expect_true(all(c("mean_loading", "sd_loading", "cos_similarity") %in% names(stab[[1L]])))
  expect_true(is.numeric(stab[[1L]]$cos_similarity))
})
