# Visualization And Interpretation Convenience Function Plan

## Purpose

This document defines a comprehensive, implementation-ready roadmap for
interpretation and visualization helpers in `rajiveplus`, with emphasis on
individual component interpretation.

Document boundary:
- This file owns API scope, sequencing, return contracts, and test planning.
- Statistical caveats and inferential warnings are authoritative in
  `StatisticalAudits.md` and should be referenced from implementation work
  rather than duplicated here.

Goals:
- Convert common post-Rajive manual analysis steps into reusable APIs.
- Provide publication-ready plots for joint and individual structure.
- Standardize statistical association workflows (clinical/phenotype linkage).
- Keep outputs tidy and composable with `ggplot2`, `dplyr`, and reporting tools.

Non-goals for initial phases:
- Rewriting decomposition internals.
- Interactive web apps as mandatory core deliverables.

## Design Principles

- Every function returns plain data (`data.frame`/list) and/or a `ggplot` object.
- Use predictable naming: `get_*` (extract), `summarize_*` (aggregate),
	`plot_*` (visual), `associate_*` (inferential), `export_*` (report output).
- Validate inputs aggressively with clear error messages.
- Preserve block/component conventions already used in the package.
- Prefer deterministic output ordering for reproducibility.
- When a planned function has methodological caveats, implement the warning by
  referencing the relevant item in `StatisticalAudits.md` instead of restating
  the full rationale here.

## Cross-File Dependencies

- `PLANS.md` decides when this roadmap becomes the active implementation lane.
- `PROGRESS.md` records execution of milestones and validation outcomes.
- `StatisticalAudits.md` defines required guardrails for:
  - `assess_stability(target = "loadings")`
  - `associate_components()`
  - `associate_scores_survival()` and split-based survival displays
  - rank-diagnostic interpretation language

## Streamlined API Surface (Recommended)

To reduce API sprawl, group similar tasks into a small set of argument-driven
functions. Keep detailed names as aliases/wrappers for discoverability.

### Core unified functions

1. `extract_components(ajive_output, blocks = NULL, source = c("rajive", "jackstraw"), what = c("scores", "loadings", "variance", "significance", "rank_diagnostics"), type = c("joint", "individual", "both"), block = NULL, component = NULL, format = c("wide", "long"), include_meta = FALSE)`
- Unifies extraction and tidy conversion (`get_*`, `tidy_*`) into one function.

2. `rank_features(ajive_output, blocks = NULL, jackstraw_result = NULL, mode = c("top_loadings", "contribution", "overlap", "significant"), type = c("joint", "individual"), block = NULL, component = NULL, n = 20, signed = TRUE, method = c("abs_loading", "variance_contrib", "jaccard", "overlap_coef"))`
- Unifies ranking/contribution/overlap/significant-feature summaries.

3. `associate_components(ajive_output, metadata, scores = NULL, mode = c("continuous", "categorical", "survival", "batch"), variable = NULL, variables = NULL, method = NULL, adjust = "BH", time_col = NULL, status_col = NULL, split = c("median", "tertile"), type = c("joint", "individual", "both"), block = NULL, component = NULL)`
- Unifies continuous/categorical/survival association wrappers.

4. `plot_components(ajive_output = NULL, blocks = NULL, jackstraw_result = NULL, assoc_results = NULL, plot_type = c("pairs", "density", "top_features", "component_heatmap", "variance", "association", "volcano", "jackstraw_summary", "stability", "ajive_diagnostic", "rank_threshold", "bound_distributions"), type = c("joint", "individual", "both"), block = NULL, component = NULL, group = NULL, style = NULL, top_n = 20, ...)`
- Unifies major plotting entry points.

5. `assess_stability(ajive_output = NULL, blocks, initial_signal_ranks, target = c("joint_rank", "loadings", "components"), method = c("bootstrap", "permutation"), B = 100, n_perm = 100, sample_frac = 0.8, num_cores = 1L)`
- Unifies stability and uncertainty diagnostics.

6. `summarize_components(ajive_output = NULL, jackstraw_result = NULL, assoc_results = NULL, summary_type = c("components", "variance", "significance", "associations", "stability", "diagnostics"), block = NULL, component = NULL, top_n = NULL)`
- Central summary table generator.

7. `export_results(x, path, format = c("csv", "tsv", "html", "md", "rds"), report = c("top_features", "component_summary", "associations", "full"), include_plots = TRUE)`
- Unifies export/report helpers.

8. `rajive_report(ajive_output, blocks, metadata = NULL, jackstraw_result = NULL, output_file = "rajive_interpretation_report.html", sections = c("overview", "variance", "features", "associations", "stability"))`
- Single-call interpretation report.

### Compatibility aliases (thin wrappers)

Retain user-friendly aliases that call unified APIs internally:
- `get_top_loadings` -> `rank_features(mode = "top_loadings")`
- `get_feature_contributions` -> `rank_features(mode = "contribution")`
- `compare_feature_sets_across_blocks` -> `rank_features(mode = "overlap")`
- `summarize_significant_vars` -> `rank_features(mode = "significant")`
- `associate_scores_continuous` -> `associate_components(mode = "continuous")`
- `associate_scores_categorical` -> `associate_components(mode = "categorical")`
- `associate_scores_survival` -> `associate_components(mode = "survival")`
- `bootstrap_joint_rank` -> `assess_stability(target = "joint_rank")`
- `bootstrap_loading_stability` -> `assess_stability(target = "loadings")`
- `plot_stability_heatmap` -> `plot_components(plot_type = "stability")`

### Streamlining rules

- New feature proposals should first attempt to extend these unified functions
	via new `mode` or `plot_type` values.
- Add a new standalone function only if argument branching becomes unclear or
	produces materially different return schemas.

## AJIVE Diagnostic Plots (Paper + mvdr alignment)

Add dedicated AJIVE rank-selection diagnostics aligned with the AJIVE paper and
the `mvdr` diagnostic style (`plot_joint_diagnostic`).

Diagnostic plot requirements:
- Visualize squared common singular values against threshold cutoffs.
- Overlay Wedin percentile cutoff and random-direction percentile cutoff.
- Use combined threshold rule: `max(wedin_cutoff, rand_cutoff)`.
- Differentiate singular values classified as joint vs non-joint.
- Optionally mark identifiability-dropped components.
- Show uncertainty bands for Wedin/random cutoff regions.

Unified API mapping:
- `plot_components(plot_type = "ajive_diagnostic")`: full composite panel.
- `plot_components(plot_type = "rank_threshold")`: focused threshold panel.
- `plot_components(plot_type = "bound_distributions")`: distributions of
	Wedin/random samples with percentile lines.

Required extraction fields:
- `extract_components(what = "rank_diagnostics")` should provide:
	- `obs_svals`
	- `wedin_samples` (optional)
	- `rand_dir_samples` (optional)
	- `wedin_cutoff` / `rand_cutoff`
	- `overall_sv_sq_threshold`
	- `joint_rank_estimate`
	- `identif_dropped` (optional)

Diagnostic acceptance criteria:
- Works with both cutoffs present, or either one missing.
- Clearly reports active cutoff in subtitle or legend.
- Returns `ggplot` and follows package theme defaults.

### Exact Function Specifications For Immediate Implementation

#### 1) `extract_components(..., what = "rank_diagnostics")`

Primary use:
- Return a normalized diagnostic payload from `ajive_output$joint_rank_sel`
  that can drive all diagnostic visualizations.

Required arguments:
- `ajive_output`: output from `Rajive()`.
- `what`: must include `"rank_diagnostics"`.

Optional arguments (supported, ignored where not relevant):
- `format`: `"wide"` (default for diagnostics) or `"long"`.
- `include_meta`: if `TRUE`, include metadata fields used in subtitles/labels.

Return schema (`format = "wide"`):
- List with class `c("rajive_diagnostics", "list")` and fields:
  - `obs_svals` numeric vector of observed singular values from common scores.
  - `obs_svals_sq` numeric vector of squared observed singular values.
  - `joint_rank_estimate` integer.
  - `overall_sv_sq_threshold` numeric scalar or `NA_real_`.
  - `wedin_samples` numeric vector or `NULL`.
  - `rand_dir_samples` numeric vector or `NULL`.
  - `wedin_cutoff` numeric scalar or `NA_real_`.
  - `rand_cutoff` numeric scalar or `NA_real_`.
  - `wedin_percentile` numeric scalar (default 5 if available).
  - `rand_percentile` numeric scalar (default 95 if available).
  - `identif_dropped` integer vector (possibly length 0).
  - `cutoff_rule` character, one of:
    - `"max(wedin, random)"`
    - `"wedin_only"`
    - `"random_only"`
    - `"none_available"`
  - `has_wedin` logical.
  - `has_random` logical.

Return schema (`format = "long"`):
- `data.frame` with columns:
  - `component_index`
  - `obs_sval`
  - `obs_sval_sq`
  - `classification` (`joint`, `nonjoint`, `dropped`)
  - `joint_rank_estimate`
  - `overall_sv_sq_threshold`
  - `wedin_cutoff`
  - `rand_cutoff`

Validation and errors:
- Error if `ajive_output$joint_rank_sel` is missing.
- Error if observed singular values are unavailable.
- Warning (not error) if both Wedin and random samples are absent.

#### 2) `plot_components(..., plot_type = "ajive_diagnostic")`

Primary use:
- Produce the full AJIVE diagnostic panel equivalent to thresholding plots in
  AJIVE references and `mvdr` style.

Required arguments:
- `ajive_output` OR precomputed diagnostics object from
  `extract_components(..., what = "rank_diagnostics")` passed via `...` as
  `diagnostics`.
- `plot_type = "ajive_diagnostic"`.

Supported `...` arguments for diagnostics:
- `diagnostics`: precomputed diagnostics payload.
- `show_wedin = TRUE`
- `show_random = TRUE`
- `show_bands = TRUE`
- `show_identif_dropped = TRUE`
- `x_scale = c("squared", "raw")`, default `"squared"`
- `wedin_percentile = 5`
- `rand_percentile = 95`
- `title = NULL`
- `subtitle = NULL`
- `theme_fn = ggplot2::theme_bw`

Plot behavior:
- Draw vertical markers for observed singular values.
- Color by classification (`joint`, `nonjoint`, `dropped`).
- Draw cutoff lines for Wedin and random bounds when available.
- Draw an emphasized line for active threshold (combined rule).
- Optionally shade percentile bands from sampled distributions.
- Add legend labels reporting cutoff values and percentiles.

Return:
- `ggplot` object.

Failure/degradation rules:
- If only one bound exists, draw only that bound and use it as active cutoff.
- If no bounds exist, draw observed values only and annotate rank as
  `"diagnostic bounds unavailable"`.
- Never fail solely because one bound source is missing.

Minimal example target calls:
- `plot_components(ajive_output = fit, plot_type = "ajive_diagnostic")`
- `plot_components(ajive_output = fit, plot_type = "rank_threshold")`
- `plot_components(ajive_output = fit, plot_type = "bound_distributions")`

## AJIVE Diagnostic Test Case Matrix

Tests cover `extract_components(what = "rank_diagnostics")` and
`plot_components(plot_type = "ajive_diagnostic" | "rank_threshold" | "bound_distributions")`.
All tests use `ajive.data.sim()` or minimal hand-constructed fixtures unless noted.

Legend — Expected column: ✓ = passes/returns, ✗ = errors with informative message, ~ = degrades gracefully with annotation.

### A. `extract_components(..., what = "rank_diagnostics")` — Return Schema Tests

| ID | Scenario | Setup | Expected |
|----|----------|-------|----------|
| A1 | Standard fit, both bounds present | `fit <- Rajive(blocks, initial_ranks)` with no override; `joint_rank_sel` contains `wedin` and `rand_dir` | Returns list with class `c("rajive_diagnostics","list")`; all 7 fields present; `joint_rank_estimate == fit$joint_rank` |
| A2 | Wide format (default) | `extract_components(fit, what="rank_diagnostics", format="wide")` | Returns data frame with one row per observed singular value; columns: `sval_index`, `obs_sval`, `wedin_cutoff`, `rand_dir_cutoff`, `overall_threshold`, `is_joint` |
| A3 | Long format | `extract_components(fit, what="rank_diagnostics", format="long")` | Returns data frame with columns `name`, `value`, `type`; `type` distinguishes `"observed"`, `"wedin_sample"`, `"rand_dir_sample"` |
| A4 | Only Wedin samples present | Set `fit$joint_rank_sel$rand_dir <- NULL` after fitting | Returns successfully; `rand_dir_cutoff = NA`; `overall_threshold = wedin_cutoff`; no error |
| A5 | Only random-direction samples present | Set `fit$joint_rank_sel$wedin <- NULL` after fitting | Returns successfully; `wedin_cutoff = NA`; `overall_threshold = rand_dir_cutoff`; no error |
| A6 | Neither bound present | Set both to NULL | Returns list with `wedin_cutoff = NA`, `rand_dir_cutoff = NA`, `overall_threshold = NA`; `joint_rank_estimate` taken directly from `$joint_rank` |
| A7 | Joint rank = 0 (all svals below threshold) | Force threshold above all observed svals | `is_joint` column is all `FALSE`; `joint_rank_estimate = 0` |
| A8 | Joint rank = full initial rank (all above) | Force threshold below all observed svals | `is_joint` column is all `TRUE`; `joint_rank_estimate` equals number of initial svals |
| A9 | Missing `joint_rank_sel` field (old object) | Remove `fit$joint_rank_sel` | Errors with message referencing `joint_rank_sel` not found; no crash |
| A10 | Single block | `blocks = list(X)`, `initial_ranks = 5` | Returns normally; `block_wedin_cutoffs` has length 1 |

### B. `plot_components(..., plot_type = "rank_threshold")` — Threshold Panel Only

| ID | Scenario | Setup | Expected |
|----|----------|-------|----------|
| B1 | Both bounds present, rand dominates | `rand_dir_cutoff > wedin_cutoff` | Random direction cutoff line drawn at `lw = 4`; Wedin line at `lw = 1`; vertical annotation at overall threshold |
| B2 | Both bounds present, Wedin dominates | `wedin_cutoff > rand_dir_cutoff` | Wedin line drawn at `lw = 4`; random direction at `lw = 1` |
| B3 | Both bounds equal | `wedin_cutoff == rand_dir_cutoff` | Both at `lw = 4` |
| B4 | Only Wedin present | `rand_dir = NULL` | Only Wedin cutoff drawn; no rand_dir line; annotation `"rand_dir bound unavailable"` present |
| B5 | Only rand present | `wedin = NULL` | Only rand_dir cutoff drawn; annotation `"Wedin bound unavailable"` present |
| B6 | No bounds at all | Both NULL | Singular values plotted; overall threshold absent; annotation `"diagnostic bounds unavailable"` present; plot still returns |
| B7 | Joint svals colored distinctly | `joint_rank = 2` | First 2 points black/solid; remaining grey/open |
| B8 | Returns ggplot | Any valid fit | `inherits(p, "ggplot")` is TRUE |

### C. `plot_components(..., plot_type = "bound_distributions")` — Distribution Histograms Only

| ID | Scenario | Setup | Expected |
|----|----------|-------|----------|
| C1 | Both distributions present | Standard fit | Two panels: Wedin bound histogram and random direction histogram; vertical line at respective cutoff (default 95th percentile) in each |
| C2 | Only Wedin distribution | `rand_dir = NULL` | Single Wedin histogram panel; no error; annotation notes missing rand_dir |
| C3 | Only rand distribution | `wedin = NULL` | Single rand_dir histogram panel; no error |
| C4 | Custom percentile cutoff | `cutoff_quantile = 0.99` | Cutoff line moves to 99th percentile; `wedin_cutoff` and `rand_dir_cutoff` in extracted diagnostics shift accordingly |
| C5 | Returns ggplot | Any valid fit | `inherits(p, "ggplot")` is TRUE |

### D. `plot_components(..., plot_type = "ajive_diagnostic")` — Full Composite Panel

| ID | Scenario | Setup | Expected |
|----|----------|-------|----------|
| D1 | Standard full composite | Standard fit with both bounds | Three-panel composite: observed svals with threshold (top), Wedin histogram (bottom-left), rand_dir histogram (bottom-right) |
| D2 | Only Wedin | `rand_dir = NULL` | Two-panel composite: svals + threshold (top), Wedin histogram (bottom); rand_dir panel absent |
| D3 | Only rand | `wedin = NULL` | Two-panel composite: svals + threshold (top), rand_dir histogram (bottom) |
| D4 | No bounds | Both NULL | Single panel: svals only with `"diagnostic bounds unavailable"` annotation |
| D5 | Returns ggplot/patchwork | Any valid fit | Result is a ggplot-compatible object; `inherits(p, "gg")` is TRUE |
| D6 | Joint rank = 0 | Forced threshold above all svals | All sval points grey; annotation shows `Joint rank: 0` |
| D7 | identif_dropped present | Use fit with identifiability drops logged | Dropped svals rendered with X marker and grey color; legend entry `"identif. dropped"` present |

### E. Error Handling And Input Validation

| ID | Scenario | Expected error message fragment |
|----|----------|---------------------------------|
| E1 | `what = "rank_diagnostics"` on jackstraw object | `"rank_diagnostics requires a rajive object"` |
| E2 | `plot_type = "ajive_diagnostic"` with no `ajive_output` | `"ajive_output must be provided"` |
| E3 | `ajive_output` is not class `"rajive"` | `"ajive_output must be of class 'rajive'"` |
| E4 | `cutoff_quantile` outside (0, 1) | `"cutoff_quantile must be in (0, 1)"` |
| E5 | `format` argument mis-spelled | `"format must be one of"` |

### F. Regression / Snapshot Tests

| ID | Scenario | Method |
|----|----------|--------|
| F1 | Canonical output for standard 2-block sim | `expect_snapshot()` on `extract_components(fit, what="rank_diagnostics", format="wide")` |
| F2 | Plot visual regression | `vdiffr::expect_doppelganger("ajive_diagnostic_standard", plot_components(fit, plot_type = "ajive_diagnostic"))` |
| F3 | Long format round-trip | `expect_equal(wide_to_long(wide_diag), long_diag)` via reshape |

## Proposed Function Catalog (Reorganized)

The catalog is grouped by user workflow from data extraction to reporting.
Each group contains high-impact functions first.

### Group 1: Data Model And Introspection Foundation

1. `get_component_info(ajive_output)`
2. `get_component_labels(ajive_output, prefix = c("J", "I"))`
3. `get_component_reconstruction_error(ajive_output, blocks, block, type)`

### Group 2: Extraction And Tidy Interoperability

4. `get_individual_scores(ajive_output, block, component = NULL)`
5. `get_individual_loadings(ajive_output, block, component = NULL)`
6. `get_joint_loadings(ajive_output, block, component = NULL)`
7. `tidy_scores(ajive_output, type = c("joint", "individual", "both"), block = NULL, sample_ids = NULL)`
8. `tidy_loadings(ajive_output, type = c("joint", "individual", "both"), block = NULL, feature_names = NULL)`
9. `tidy_variance_explained(ajive_output, blocks)`
10. `tidy_jackstraw(jackstraw_result, block = NULL, component = NULL)`

### Group 3: Individual Component Interpretation Core

11. `get_individual_top_features(ajive_output, block, component, n = 20, direction = c("both", "positive", "negative"), abs_rank = TRUE)`
12. `summarize_individual_top_features(ajive_output, block, n = 20, method = c("frequency", "mean_abs_loading"))`
13. `get_individual_extreme_samples(ajive_output, block, component, method = c("zscore", "quantile"), threshold = 2.5, q = 0.95)`
14. `summarize_individual_variance(ajive_output, blocks, block = NULL, cumulative = TRUE)`
15. `get_individual_feature_sets(ajive_output, block, top_n = 50, overlap = c("jaccard", "overlap_coef"))`
16. `rank_individual_features_by_stability(ajive_output, blocks, block, component, B = 50, sample_frac = 0.8)`

### Group 4: Feature Contribution And Ranking Utilities

17. `get_top_loadings(ajive_output, block, component, type = c("joint", "individual"), n = 20, signed = TRUE)`
18. `get_feature_contributions(ajive_output, blocks, block, component, type = c("joint", "individual"), method = c("abs_loading", "variance_contrib"))`
19. `compare_feature_sets_across_blocks(ajive_output, component, type = c("joint", "individual"), top_n = 50, overlap = c("jaccard", "overlap_coef"))`

### Group 5: Metadata Association And Inference

20. `associate_scores_continuous(scores, metadata, variable, method = c("pearson", "spearman"), adjust = "BH")`
21. `associate_scores_categorical(scores, metadata, variable, method = c("anova", "kruskal"), adjust = "BH")`
22. `associate_scores_survival(scores, metadata, time_col, status_col, model = c("coxph", "km_split"), split = c("median", "tertile"))`
23. `associate_components(ajive_output, metadata, types = c("joint", "individual"), blocks = NULL, tests = c("numeric", "factor"), adjust = "BH")`
24. `summarize_associations(assoc_results, top_n = 20, by = c("p_value", "effect_size"))`

### Group 6: Visualization Layer (Unified + Specialized)

Unified visualization entry points:
25. `autoplot.rajive(object, plot_type = c("variance", "scores", "pairs", "top_features", "component_heatmap"), ...)`
26. `autoplot.jackstraw_rajive(object, plot_type = c("hist", "scatter", "heatmap", "summary", "volcano"), ...)`

Core plot primitives:
27. `plot_component_pairs(ajive_output, type = c("joint", "individual"), block = NULL, group = NULL)`
28. `plot_individual_score_density(ajive_output, block, component, group = NULL, ridge = FALSE)`
29. `plot_individual_top_features(ajive_output, block, component, n = 20, direction = "both", style = c("bar", "lollipop"))`
30. `plot_individual_feature_heatmap(block_matrix, ajive_output, block, component, top_n = 40, order_by_score = TRUE)`
31. `plot_joint_scores(ajive_output, comp_x = 1, comp_y = 2, group = NULL)`
32. `plot_individual_scores(ajive_output, block, comp_x = 1, comp_y = 2, group = NULL, ellipse = FALSE)`
33. `plot_variance_waterfall(ajive_output, blocks, block = NULL)`
34. `plot_component_contribution_map(ajive_output, blocks, normalize = c("block", "global"))`
35. `plot_component_correlation(ajive_output, metadata = NULL, method = c("pearson", "spearman"))`
36. `plot_individual_feature_overlap(ajive_output, block, top_n = 50, method = c("heatmap", "network"))`
37. `plot_individual_sample_outliers(ajive_output, block, component, threshold = 2.5, label = TRUE)`

Association and significance visuals:
38. `plot_component_associations(assoc_results, style = c("forest", "heatmap", "lollipop"))`
39. `plot_survival_by_component(scores, metadata, time_col, status_col, component, split = "median")`
40. `plot_component_effect_size(assoc_results, metric = c("r", "eta2", "HR"))`
41. `plot_jackstraw_volcano(jackstraw_result, block, component)`
42. `plot_jackstraw_component_summary(jackstraw_result, block = NULL)`
43. `plot_significant_feature_upset(jackstraw_result, block, top_n = 30)`
44. `plot_significance_vs_loading(jackstraw_result, ajive_output, block, component)`

### Group 7: Jackstraw Summaries And Convenience Accessors

45. `summarize_significant_vars(jackstraw_result, block = NULL, component = NULL, top_n = NULL)`

### Group 8: Stability And Uncertainty Diagnostics

46. `bootstrap_joint_rank(blocks, initial_signal_ranks, B = 100, sample_frac = 0.8, num_cores = 1L)`
47. `bootstrap_loading_stability(ajive_output, blocks, block, component, B = 50, sample_frac = 0.8, num_cores = 1L)`
48. `bootstrap_component_stability(blocks, initial_signal_ranks, B = 50, sample_frac = 0.8, num_cores = 1L)`
49. `summarize_component_stability(stability_result)`
50. `plot_stability_heatmap(stability_result)`
51. `plot_component_stability(stability_result, style = c("heatmap", "line", "interval"))`
52. `permute_component_significance(ajive_output, blocks, block, component, n_perm = 100)`

### Group 9: Reporting And Export

53. `export_top_features_csv(df, path)`
54. `export_component_summary_html(summary_obj, path)`
55. `export_component_table(df, path, format = c("csv", "tsv"))`
56. `export_top_features_report(ajive_output, path, format = c("html", "md"))`
57. `export_association_report(assoc_results, path, include_plots = TRUE)`
58. `rajive_report(ajive_output, blocks, metadata = NULL, output_file = "rajive_interpretation_report.html")`

### Group 10: S3 Data-Prep Methods

59. `fortify.rajive(model, data, ...)`
60. `fortify.jackstraw_rajive(model, data, ...)`

## Dependency Plan

Required (existing or lightweight):
- `ggplot2`
- `stats`

Suggested additions:
- `ggrepel` (labels)
- `pheatmap` or `ComplexHeatmap` (heatmaps)
- `survival`, `survminer` (survival association)
- `reshape2` replacement should move to `tidyr` where possible
- `dplyr` optional for internal clarity (can avoid to reduce hard deps)

Policy:
- Put plotting/report extras in `Suggests` with graceful `requireNamespace` checks.
- Keep core extract/associate helpers in strict minimal dependency footprint.

## Actionable Implementation Roadmap

### Phase 0: Spec And Infrastructure (0.5 to 1 day)

Deliverables:
- Finalize unified function signatures and allowed enum values for
  `what`, `mode`, `plot_type`, `target`, and `summary_type`.
- Freeze diagnostic enums and fields for AJIVE diagnostics:
	`rank_diagnostics`, `ajive_diagnostic`, `rank_threshold`,
	`bound_distributions`.
- Add shared validators:
	- block/component index checks
	- metadata alignment checks
	- score/loading availability checks
- Add internal helper for consistent long-table schemas and return contracts.

Acceptance criteria:
- Naming and argument conventions frozen for Phase 1.
- Internal helper unit tests pass.

### Phase 1: Individual Interpretation MVP (1 to 2 days)

Implement first:
- `extract_components` (scores/loadings/variance; long format)
- `rank_features` (top_loadings, contribution)
- `plot_components` (`top_features`, `pairs`, `density`, `component_heatmap`)
- `extract_components` (`what = "rank_diagnostics"`)
- `plot_components` (`rank_threshold`)
- thin aliases: `get_top_loadings`, `get_individual_top_features`,
  `summarize_individual_variance`

Tests:
- deterministic ranking on simulated data
- edge cases: rank-1 components, missing feature names, NA handling
- plot object classes and minimal structural checks

Docs:
- roxygen examples for each function
- one mini article section in vignette showing end-to-end usage

Acceptance criteria:
- all functions exported and documented
- tests green
- examples run in `donttest` without errors

### Phase 2: Association Layer (1 to 2 days)

Implement:
- `associate_components` (`continuous`, `categorical`, `batch`)
- `summarize_components` (`associations`)
- `plot_components` (`association`)
- thin aliases: `associate_scores_continuous`,
  `associate_scores_categorical`

Tests:
- known synthetic associations recovered
- p-value adjustment behavior verified
- metadata/sample mismatch throws actionable error

Acceptance criteria:
- tidy outputs with stable columns
- at least one vignette section using clinical metadata

### Phase 3: Jackstraw Integration Enhancements (1 day)

Implement:
- `extract_components` (`source = "jackstraw"`, `what = "significance"`)
- `rank_features` (`mode = "significant"`)
- `plot_components` (`volcano`, `jackstraw_summary`, `ajive_diagnostic`,
  `bound_distributions`)
- thin alias: `summarize_significant_vars`

Acceptance criteria:
- compatibility with current `jackstraw_rajive` object structure
- plot functions return `ggplot` objects consistently

### Phase 4: Stability Toolkit (1 to 2 days)

Implement:
- `assess_stability` (`target = "joint_rank"`, `target = "loadings"`)
- `summarize_components` (`summary_type = "stability"`)
- `plot_components` (`plot_type = "stability"`)
- thin aliases: `bootstrap_joint_rank`, `bootstrap_loading_stability`,
  `plot_stability_heatmap`

Acceptance criteria:
- runtime-configurable bootstrap size
- reproducibility with fixed seed
- clear warning for expensive settings

### Phase 5: Reporting And S3 Ergonomics (1 day)

Implement:
- `export_results`
- `rajive_report`
- `autoplot.rajive` and `autoplot.jackstraw_rajive` as wrappers around
	`plot_components`

Acceptance criteria:
- single-call summary output possible for common workflows
- no hard dependency on heavy optional packages

## Plan Audit

### Audit Outcome

Status: pass with targeted refinements.

What is strong:
- End-to-end workflow coverage from extraction to reporting.
- Clear balance of interpretation, statistics, and visualization.
- Explicit support for both joint and individual component workflows.
- Practical phased rollout with testing expectations.

Issues found and fixes applied:
1. Naming mismatch risk between requested APIs and planned APIs.
- Fix: added explicit requested names as first-class functions:
	`get_top_loadings`, `get_feature_contributions`,
	`compare_feature_sets_across_blocks`, `associate_scores_continuous`,
	`associate_scores_categorical`, `associate_scores_survival`,
	`summarize_significant_vars`, `bootstrap_joint_rank`,
	`bootstrap_loading_stability`, `plot_stability_heatmap`,
	`export_top_features_csv`, `export_component_summary_html`, `rajive_report`.

2. Discoverability risk due to broad function surface.
- Fix: reorganized into 10 logical groups following user workflow order.

3. Potential over-investment in low-usage advanced features early.
- Fix: tightened priority tiers to emphasize high-impact interpretation functions
	before advanced stability/reporting features.

4. Redundancy risk between unified plotting and specific plotting functions.
- Fix: defined unified entry points (`autoplot.*`) and retained specific
	plot primitives as composable internals.

Residual risks:
- API size may still feel large for first-time users.
- Stability functions may be compute-intensive on large datasets.
- Survival helpers require optional dependencies and robust metadata checks.

Mitigations:
- ship a compact Phase 1 interface first and expose advanced functions later
- add strict input validators and informative warnings
- keep heavy dependencies in Suggests and fail gracefully

## File-Level Execution Plan

Suggested file split:
- `R/interpretation_individual.R`
- `R/interpretation_association.R`
- `R/interpretation_tidy.R`
- `R/visualization_individual.R`
- `R/visualization_jackstraw_ext.R`
- `R/stability.R`
- `R/reporting.R`

Tests:
- `tests/testthat/test-interpretation-individual.R`
- `tests/testthat/test-interpretation-association.R`
- `tests/testthat/test-interpretation-tidy.R`
- `tests/testthat/test-visualization-individual.R`
- `tests/testthat/test-visualization-jackstraw-ext.R`
- `tests/testthat/test-stability.R`
- `tests/testthat/test-reporting.R`

Documentation:
- add reference index sections in `_pkgdown.yml`:
	- Interpretation
	- Visualization
	- Association
	- Stability

## Milestone Backlog With Priority

Priority P0: highest impact, fastest user value
- `extract_components`
- `rank_features`
- `plot_components`
- `autoplot.rajive`
- compatibility aliases for top loadings and individual summaries
- `plot_components(plot_type = "rank_threshold")`

Priority P1: interpretation-to-biology bridge
- `associate_components`
- `summarize_components` (associations and significance)
- `plot_components` (association and volcano)
- compatibility aliases for continuous/categorical/survival helpers
- `plot_components(plot_type = "ajive_diagnostic")`
- `plot_components(plot_type = "bound_distributions")`

Priority P2: credibility and robustness
- `assess_stability`
- `summarize_components` (stability)
- `plot_components` (stability)

Priority P3: reporting polish and advanced extras
- `export_results`
- `rajive_report`
- advanced aliases and wrappers
- `fortify.rajive`
- `fortify.jackstraw_rajive`

Recommended first 10-function implementation bundle:
1. `extract_components`
2. `rank_features`
3. `plot_components`
4. `associate_components`
5. `summarize_components`
6. `assess_stability`
7. `export_results`
8. `rajive_report`
9. `autoplot.rajive`
10. `autoplot.jackstraw_rajive`

## Definition Of Done (Per Function)

Each function is considered complete only if all are true:
- exported and documented with examples
- covered by unit tests including failure paths
- integrated into at least one vignette/example workflow
- included in NEWS (if release-facing)
- appears in pkgdown reference site with category placement

## Immediate Next Action

Start Phase 1 with the individual interpretation MVP:

1. Implement `extract_components` and `rank_features` with strict return schemas.
2. Implement `plot_components` for `pairs`, `density`, `top_features`, and `component_heatmap`.
2.1 Add AJIVE diagnostics: `rank_threshold` first, then `ajive_diagnostic` and `bound_distributions`.
3. Add compatibility aliases (`get_top_loadings`, `get_individual_top_features`, `tidy_scores`, `tidy_loadings`) as thin wrappers.
4. Add tests and minimal vignette examples before proceeding to association integration.

## Requested Suggestions Coverage Checklist

This section maps requested suggestions to the streamlined unified API.

Unified plotting API:
- Requested: autoplot methods with plot_type argument
- Planned: `autoplot.rajive` and `autoplot.jackstraw_rajive` (wrappers over `plot_components`)
- Requested: score_pairs
- Planned: `plot_components(plot_type = "pairs")`
- Requested: score_density by group
- Planned: `plot_components(plot_type = "density", group = ...)`
- Requested: loading_barplot
- Planned: `plot_components(plot_type = "top_features", style = "bar")`
- Requested: component_heatmap for sample scores
- Planned: `plot_components(plot_type = "component_heatmap")`
- Requested: AJIVE diagnostic threshold plot (Wedin + random direction)
- Planned: `plot_components(plot_type = "ajive_diagnostic")`,
  `plot_components(plot_type = "rank_threshold")`, and
  `plot_components(plot_type = "bound_distributions")`

Feature contribution and ranking helpers:
- Requested: get_top_loadings(block, component, type, n, signed)
- Planned: alias `get_top_loadings` -> `rank_features(mode = "top_loadings")`
- Requested: get_feature_contributions(method = abs_loading or variance_contrib)
- Planned: alias `get_feature_contributions` -> `rank_features(mode = "contribution")`
- Requested: compare_feature_sets_across_blocks(component)
- Planned: alias `compare_feature_sets_across_blocks` -> `rank_features(mode = "overlap")`

Metadata association module:
- Requested: associate_scores_continuous
- Planned: alias `associate_scores_continuous` -> `associate_components(mode = "continuous")`
- Requested: associate_scores_categorical
- Planned: alias `associate_scores_categorical` -> `associate_components(mode = "categorical")`
- Requested: associate_scores_survival
- Planned: alias `associate_scores_survival` -> `associate_components(mode = "survival")`
- Requested: plot_component_associations
- Planned: `plot_components(plot_type = "association")`

Jackstraw downstream summaries:
- Requested: as_tibble_jackstraw
- Planned: `extract_components(source = "jackstraw", what = "significance", format = "long")`
- Requested: summarize_significant_vars by block/component
- Planned: alias `summarize_significant_vars` -> `rank_features(mode = "significant")`
- Requested: plot_significance_upset
- Planned: `plot_components(plot_type = "jackstraw_summary", style = "upset")`

Stability and uncertainty diagnostics:
- Requested: bootstrap_joint_rank
- Planned: alias `bootstrap_joint_rank` -> `assess_stability(target = "joint_rank")`
- Requested: bootstrap_loading_stability
- Planned: alias `bootstrap_loading_stability` -> `assess_stability(target = "loadings")`
- Requested: plot_stability_heatmap
- Planned: alias `plot_stability_heatmap` -> `plot_components(plot_type = "stability", style = "heatmap")`

Report-ready exports:
- Requested: export_top_features_csv
- Planned: alias `export_top_features_csv` -> `export_results(report = "top_features", format = "csv")`
- Requested: export_component_summary_html
- Planned: alias `export_component_summary_html` -> `export_results(report = "component_summary", format = "html")`
- Requested: rajive_report
- Planned: `rajive_report(...)`

Interoperability helpers:
- Requested: tidy_scores
- Planned: alias `tidy_scores` -> `extract_components(what = "scores", format = "long")`
- Requested: tidy_loadings
- Planned: alias `tidy_loadings` -> `extract_components(what = "loadings", format = "long")`
- Requested: tidy_variance_explained
- Planned: alias `tidy_variance_explained` -> `extract_components(what = "variance", format = "long")`

Coverage status:
- All requested capabilities are now covered in the streamlined surface,
  with explicit alias paths where helpful for readability.

