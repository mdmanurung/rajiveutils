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
                       "identif_dropped", "identifiability_norm", "cutoff_rule",
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

# ---------------------------------------------------------------------------
# H. rank_features — return schema and correctness tests
# ---------------------------------------------------------------------------

# Shared jackstraw fixture (lazy, computed once)
.viz_jsr <- local({
  set.seed(99L)
  jackstraw_rajive(.viz_fit, ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                                            n = 40L, pks = c(30L, 20L),
                                            dist.type = 1)$sim_data,
                   n_null = 5L)
})

test_that("H1: top_loadings returns correct columns and respects n", {
  fit <- .viz_fit
  out <- rank_features(fit, mode = "top_loadings", type = "joint", n = 5L)

  expect_true(is.data.frame(out))
  expect_true(all(c("block", "component", "type", "feature_index",
                    "feature_name", "loading", "abs_loading", "rank") %in% names(out)))
  # n = 5 per block x component combination
  expect_true(all(out$rank <= 5L))
  # ranking is by abs_loading descending within each block x component
  for (bc in unique(paste(out$block, out$component))) {
    sub <- out[paste(out$block, out$component) == bc, ]
    expect_true(all(diff(sub$abs_loading) <= 0))
  }
})

test_that("H2: top_loadings signed=TRUE preserves sign; signed=FALSE gives abs values", {
  fit     <- .viz_fit
  out_s   <- rank_features(fit, mode = "top_loadings", type = "individual",
                            block = 1L, component = 1L, n = 10L, signed = TRUE)
  out_u   <- rank_features(fit, mode = "top_loadings", type = "individual",
                            block = 1L, component = 1L, n = 10L, signed = FALSE)

  expect_true(all(out_u$loading >= 0))
  expect_equal(abs(out_s$loading), out_s$abs_loading)
  expect_equal(out_u$loading, out_s$abs_loading)
})

test_that("H3: contribution mode adds contribution column, contributions are positive", {
  fit <- .viz_fit
  out <- rank_features(fit, mode = "contribution", type = "joint",
                       block = 1L, component = 1L, n = 10L,
                       method = "variance_contrib")

  expect_true("contribution" %in% names(out))
  expect_true(all(out$contribution >= 0))
  expect_true(all(out$contribution <= 1))
  # top feature has highest contribution
  expect_true(all(diff(out$contribution) <= 0))
})

test_that("H4: contribution abs_loading method contributions are positive and rank-ordered", {
  fit <- .viz_fit
  out <- rank_features(fit, mode = "contribution", type = "individual",
                       block = 1L, n = 5L, method = "abs_loading")

  expect_true(all(out$contribution >= 0))
  for (bc in unique(paste(out$block, out$component))) {
    sub <- out[paste(out$block, out$component) == bc, ]
    expect_true(all(diff(sub$contribution) <= 0))
  }
})

test_that("H5: overlap mode returns correct columns and scores in [0, 1]", {
  fit <- .viz_fit
  out <- rank_features(fit, mode = "overlap", type = "joint",
                       component = 1L, n = 10L, method = "jaccard")

  expect_true(is.data.frame(out))
  expect_true(all(c("component", "type", "block_i", "block_j",
                    "top_n", "method", "n_intersect", "overlap_score") %in% names(out)))
  expect_true(all(out$overlap_score >= 0 & out$overlap_score <= 1, na.rm = TRUE))
  expect_true(all(out$block_i < out$block_j))  # upper triangle
})

test_that("H6: overlap overlap_coef method gives valid scores", {
  fit <- .viz_fit
  out <- rank_features(fit, mode = "overlap", type = "individual",
                       n = 5L, method = "overlap_coef")

  expect_true(all(out$overlap_score >= 0 & out$overlap_score <= 1, na.rm = TRUE))
})

test_that("H7: significant mode returns correct columns for jackstraw result", {
  jsr <- .viz_jsr
  out <- rank_features(jackstraw_result = jsr, mode = "significant")

  expect_true(is.data.frame(out))
  expect_true(all(c("block", "component", "feature_index",
                    "p_value", "p_adj", "significant") %in% names(out)))
  expect_true(all(out$p_value >= 0 & out$p_value <= 1))
  expect_true(is.logical(out$significant))
})

test_that("H8: significant mode block/component filtering works", {
  jsr    <- .viz_jsr
  out_all  <- rank_features(jackstraw_result = jsr, mode = "significant")
  out_b1   <- rank_features(jackstraw_result = jsr, mode = "significant",
                             block = 1L)
  out_b2   <- rank_features(jackstraw_result = jsr, mode = "significant",
                             block = 2L)

  expect_true(all(out_b1$block == 1L))
  expect_true(all(out_b2$block == 2L))
  expect_equal(nrow(out_all), nrow(out_b1) + nrow(out_b2))
})

test_that("H9: n larger than feature count returns all features", {
  fit <- .viz_fit
  # block 2 has 20 features; requesting 100 should return all 20
  out <- rank_features(fit, mode = "top_loadings", type = "joint",
                       block = 2L, component = 1L, n = 100L)
  # joint loadings for block 2 have 20 rows
  expect_equal(nrow(out), 20L)
})

test_that("H10: invalid method for mode raises error", {
  fit <- .viz_fit
  expect_error(
    rank_features(fit, mode = "top_loadings", method = "jaccard"),
    regexp = "not valid for mode"
  )
  expect_error(
    rank_features(fit, mode = "overlap", method = "abs_loading"),
    regexp = "not valid for mode"
  )
})

test_that("H11: significant mode errors without jackstraw_result", {
  expect_error(
    rank_features(mode = "significant"),
    regexp = "jackstraw_result.*must be provided"
  )
})

test_that("H12: aliases delegate correctly to rank_features", {
  fit <- .viz_fit

  tl  <- get_top_loadings(fit, block = 1L, component = 1L, n = 5L)
  fc  <- get_feature_contributions(fit, block = 1L, component = 1L)
  ov  <- compare_feature_sets_across_blocks(fit, component = 1L, n = 10L)
  sv  <- summarize_significant_vars(.viz_jsr, block = 1L)

  expect_true(is.data.frame(tl))
  expect_true("rank" %in% names(tl))
  expect_true("contribution" %in% names(fc))
  expect_true("overlap_score" %in% names(ov))
  expect_true("p_value" %in% names(sv))
  expect_true(all(sv$block == 1L))
})

test_that("H13: summarize_significant_vars top_n argument limits rows per group", {
  jsr <- .viz_jsr
  sv  <- summarize_significant_vars(jsr, top_n = 3L)

  counts <- table(paste(sv$block, sv$component))
  expect_true(all(counts <= 3L))
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


# ---------------------------------------------------------------------------
# I. summarize_components — summary table tests
# ---------------------------------------------------------------------------

test_that("I1: summarize_components('components') returns expected schema", {
  out <- summarize_components(ajive_output = .viz_fit, summary_type = "components")
  expect_true(is.data.frame(out))
  expect_true(all(c("block", "type", "component", "block_rank",
                    "joint_rank", "n_features") %in% names(out)))
  expect_true(all(out$type %in% c("joint", "individual")))
})

test_that("I2: summarize_components('diagnostics') mirrors extract long diagnostics", {
  out <- summarize_components(ajive_output = .viz_fit, summary_type = "diagnostics")
  expect_true(is.data.frame(out))
  expect_true(all(c("component_index", "obs_sval", "classification") %in% names(out)))
})

test_that("I3: summarize_components('significance') returns per block-component aggregates", {
  out <- summarize_components(jackstraw_result = .viz_jsr, summary_type = "significance")
  expect_true(is.data.frame(out))
  expect_true(all(c("block", "component", "n_features",
                    "n_significant", "prop_significant") %in% names(out)))
  expect_true(all(out$prop_significant >= 0 & out$prop_significant <= 1))
})

test_that("I4: summarize_components('associations') aggregates variable-component rows", {
  n <- nrow(.viz_fit$joint_scores)
  meta <- data.frame(
    x = rnorm(n),
    y = sample(c("A", "B"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  assoc <- suppressMessages(
    associate_components(.viz_fit, meta, variables = c("x", "y"), mode = "categorical")
  )

  out <- summarize_components(assoc_results = assoc, summary_type = "associations")
  expect_true(is.data.frame(out))
  expect_true(all(c("variable", "component", "n_tests", "min_p_value", "min_p_adj") %in% names(out)))
})

test_that("I5: summarize_components('stability') handles rank-table output", {
  stab <- assess_stability(
    ajive_output = .viz_fit,
    blocks = ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                            n = 40L, pks = c(30L, 20L),
                            dist.type = 1)$sim_data,
    initial_signal_ranks = c(4L, 3L),
    target = "joint_rank",
    B = 3L,
    sample_frac = 0.8
  )

  out <- summarize_components(stability_result = stab, summary_type = "stability")
  expect_true(is.data.frame(out))
  expect_true(all(c("joint_rank", "frequency", "proportion") %in% names(out)))
})

test_that("I6: summarize_components validates required inputs", {
  expect_error(
    summarize_components(summary_type = "components"),
    regexp = "ajive_output"
  )
  expect_error(
    summarize_components(summary_type = "significance"),
    regexp = "jackstraw_result"
  )
  expect_error(
    summarize_components(summary_type = "associations"),
    regexp = "assoc_results"
  )
})

# ---------------------------------------------------------------------------
# J. Remaining non-benchmark modules
# ---------------------------------------------------------------------------

test_that("J1: extract_components supports scores/loadings/significance/variance", {
  Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                      n = 40L, pks = c(30L, 20L), dist.type = 1)
  fit <- .viz_fit

  s_long <- extract_components(fit, what = "scores", format = "long")
  l_long <- extract_components(fit, what = "loadings", format = "long")
  v_tbl  <- extract_components(fit, blocks = Y$sim_data, what = "variance")
  sig    <- extract_components(what = "significance", source = "jackstraw",
                               jackstraw_result = .viz_jsr)

  expect_true(all(c("score", "component") %in% names(s_long)))
  expect_true(all(c("loading", "feature") %in% names(l_long)))
  expect_true(is.data.frame(v_tbl))
  expect_true(all(c("p_value", "significant") %in% names(sig)))
})

test_that("J1b: plot_components supports loading confidence intervals", {
  fx <- make_small_rajive_fixture(seed = 8101L)
  set.seed(8102L)
  ci <- rajive_ci(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    target = "loadings",
    method = "percentile",
    B = 3L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  p <- plot_components(fx$fit, plot_type = "loading_ci", ci = ci, top_n = 5L)
  expect_true(inherits(p, "ggplot"))
})

test_that("J2: plot_components supports non-diagnostic plot types", {
  Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                      n = 40L, pks = c(30L, 20L), dist.type = 1)
  assoc <- suppressMessages(
    associate_components(.viz_fit,
                         data.frame(g = sample(c("A", "B"), 40L, TRUE), stringsAsFactors = FALSE),
                         variable = "g", mode = "categorical")
  )
  stab <- assess_stability(ajive_output = .viz_fit,
                           blocks = Y$sim_data,
                           initial_signal_ranks = c(4L, 3L),
                           target = "joint_rank",
                           B = 3L)

  expect_true(inherits(plot_components(.viz_fit, plot_type = "pairs"), "ggplot"))
  expect_true(inherits(plot_components(.viz_fit, plot_type = "density"), "ggplot"))
  expect_true(inherits(plot_components(.viz_fit, plot_type = "top_features"), "ggplot"))
  expect_true(inherits(plot_components(.viz_fit, plot_type = "component_heatmap"), "ggplot"))
  expect_true(inherits(plot_components(.viz_fit, blocks = Y$sim_data, plot_type = "variance"), "ggplot"))
  expect_true(inherits(plot_components(assoc_results = assoc, plot_type = "association"), "ggplot"))
  expect_true(inherits(plot_components(jackstraw_result = .viz_jsr, plot_type = "volcano"), "ggplot"))
  expect_true(inherits(plot_components(jackstraw_result = .viz_jsr, plot_type = "jackstraw_summary"), "ggplot"))
  expect_true(inherits(plot_components(stability_result = stab, plot_type = "stability"), "ggplot"))
})

test_that("J2b: association plots support all-component identities and uncertainty", {
  fx <- make_small_rajive_fixture(seed = 8201L)
  metadata <- data.frame(marker = fx$fit$joint_scores[, 1L] + rnorm(12, sd = 0.1))

  set.seed(8202L)
  assoc_all <- suppressMessages(
    associate_all_components(
      fx$fit, metadata, variable = "marker",
      mode = "continuous",
      include = "both",
      blocks_to_include = 1L,
      propagate_uncertainty = "bootstrap",
      blocks = fx$blocks,
      initial_signal_ranks = fx$initial_signal_ranks,
      B = 2L,
      n_wedin_samples = NA,
      n_rand_dir_samples = NA,
      joint_rank = 1L
    )
  )

  p_heat <- plot_components(assoc_results = assoc_all, plot_type = "association")
  p_forest <- plot_components(assoc_results = assoc_all, plot_type = "association",
                              style = "forest")
  p_unc <- plot_components(assoc_results = assoc_all, plot_type = "association",
                           style = "uncertainty")

  expect_true(inherits(p_heat, "ggplot"))
  expect_true(inherits(p_forest, "ggplot"))
  expect_true(inherits(p_unc, "ggplot"))
  expect_true(any(grepl("individual", p_heat$data$target)))
})

test_that("J3: export_results writes files", {
  td <- tempdir()
  p1 <- file.path(td, "exp_results.csv")
  p2 <- file.path(td, "exp_results.rds")
  p3 <- file.path(td, "exp_results.md")

  export_results(data.frame(a = 1:3), p1, format = "csv")
  export_results(list(a = 1), p2, format = "rds")
  export_results(data.frame(a = 1:3), p3, format = "md")

  expect_true(file.exists(p1))
  expect_true(file.exists(p2))
  expect_true(file.exists(p3))
})

test_that("J4: rajive_report writes output file", {
  Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                      n = 40L, pks = c(30L, 20L), dist.type = 1)
  out <- file.path(tempdir(), "rajive_report_test.md")
  rajive_report(.viz_fit, Y$sim_data, output_file = out,
                sections = c("overview", "variance", "features"))
  expect_true(file.exists(out))
})

test_that("J5: alias wrappers run and return expected object types", {
  meta <- data.frame(x = rnorm(40L), g = sample(c("A", "B"), 40L, TRUE),
                     time = rexp(40L), status = sample(0:1, 40L, TRUE),
                     stringsAsFactors = FALSE)
  Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                      n = 40L, pks = c(30L, 20L), dist.type = 1)

  ac <- suppressMessages(associate_scores_continuous(.viz_fit, meta, variable = "x"))
  ag <- suppressMessages(associate_scores_categorical(.viz_fit, meta, variable = "g"))
  as <- suppressMessages(associate_scores_survival(.viz_fit, meta,
                                                   time_col = "time", status_col = "status"))
  bj <- bootstrap_joint_rank(.viz_fit, Y$sim_data, c(4L, 3L), B = 3L)
  bl <- bootstrap_loading_stability(.viz_fit, Y$sim_data, c(4L, 3L), B = 3L)

  expect_true(is.data.frame(ac))
  expect_true(is.data.frame(ag))
  expect_true(is.data.frame(as))
  expect_true(is.list(bj))
  expect_true(is.list(bl))
})

test_that("J6: assess_stability target='components' returns summary table", {
  Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                      n = 40L, pks = c(30L, 20L), dist.type = 1)
  out <- assess_stability(ajive_output = .viz_fit,
                          blocks = Y$sim_data,
                          initial_signal_ranks = c(4L, 3L),
                          target = "components",
                          B = 3L)
  expect_true(is.data.frame(out))
  expect_true(all(c("component", "mean_correlation", "sd_correlation") %in% names(out)))
})

# ---------------------------------------------------------------------------
# K. autoplot / fortify S3 methods
# ---------------------------------------------------------------------------

test_that("K1: autoplot.rajive returns ggplot for pairs/density/top_features", {
  expect_true(inherits(ggplot2::autoplot(.viz_fit, plot_type = "pairs"),        "ggplot"))
  expect_true(inherits(ggplot2::autoplot(.viz_fit, plot_type = "density"),      "ggplot"))
  expect_true(inherits(ggplot2::autoplot(.viz_fit, plot_type = "top_features"), "ggplot"))
})

test_that("K1b: autoplot.rajive supports loading confidence intervals", {
  fx <- make_small_rajive_fixture(seed = 8201L)
  set.seed(8202L)
  ci <- rajive_ci(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    target = "loadings",
    method = "percentile",
    B = 3L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  expect_true(inherits(ggplot2::autoplot(fx$fit, plot_type = "loading_ci",
                                         ci = ci, top_n = 5L), "ggplot"))
})

test_that("K2: autoplot.rajive supports diagnostic plot_types", {
  expect_true(inherits(ggplot2::autoplot(.viz_fit, plot_type = "rank_threshold"),      "ggplot"))
  expect_true(inherits(ggplot2::autoplot(.viz_fit, plot_type = "bound_distributions"), "ggplot"))
  expect_true(inherits(ggplot2::autoplot(.viz_fit, plot_type = "ajive_diagnostic"),    "gg"))
})

test_that("K3: autoplot.jackstraw_rajive returns ggplot for jackstraw_summary", {
  expect_true(inherits(ggplot2::autoplot(.viz_jsr, plot_type = "jackstraw_summary"), "ggplot"))
})

test_that("K4: autoplot.jackstraw_rajive returns ggplot for volcano", {
  expect_true(inherits(ggplot2::autoplot(.viz_jsr, plot_type = "volcano",
                                         block = 1L, component = 1L), "ggplot"))
})

test_that("K5: fortify.rajive returns long data.frame with score column", {
  df <- ggplot2::fortify(.viz_fit, what = "scores")
  expect_true(is.data.frame(df))
  expect_true("score" %in% names(df))
  expect_true("component" %in% names(df))
})

test_that("K6: fortify.rajive supports loadings", {
  df <- ggplot2::fortify(.viz_fit, what = "loadings")
  expect_true(is.data.frame(df))
  expect_true("loading" %in% names(df))
})

test_that("K7: fortify.jackstraw_rajive returns significance data.frame", {
  df <- ggplot2::fortify(.viz_jsr)
  expect_true(is.data.frame(df))
  expect_true(all(c("p_value", "significant") %in% names(df)))
})

test_that("K8: autoplot errors informatively for wrong object class", {
  expect_error(rajiveplus:::autoplot.rajive(list()), regexp = "class")
})
