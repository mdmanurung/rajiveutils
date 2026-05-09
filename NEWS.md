# rajiveplus (development)

## New features

- `jackstraw_rajive()` gains three new arguments for posterior inclusion
  probabilities: `pip` (default `FALSE`), `pip_pi0`, and `pip_group`
  (`"component"`, `"block_component"`, or `"pooled"`).  When `pip = TRUE`,
  each block/component result list includes a `$pip` vector of per-feature
  PIPs computed as `1 - qvalue::lfdr(p_values)`.  Requires the
  **qvalue** Bioconductor package (Suggests).
  Reference: Chung NC (2020). *Bioinformatics*, 36(10):3107–3114.

## Changes in defaults

- `jackstraw_rajive()`: default `correction` changed from `"BH"` to `"BY"`
  (Benjamini--Yekutieli).  BY controls FDR under arbitrary dependency among
  tests; BH requires positive regression dependency (PRDS).  Omics feature
  blocks (with local correlation structure) can violate PRDS, making BH
  anti-conservative.  Existing code that passed `correction = "BH"` explicitly
  is unaffected.

- Added `n_perm_samples` argument to `Rajive()` (default `NA`, no behaviour
  change). When non-NA, a non-parametric permutation-based joint-rank
  threshold is computed by independently permuting the row order of each
  block's signal scores and sampling the leading squared singular value of the
  stacked matrix. This replaces the random direction bound and provides a
  data-driven null distribution that adapts to the actual noise in each block,
  making it particularly useful for heavy-tailed or non-Gaussian data.
  Diagnostic results are stored in `fit$joint_rank_sel$perm` and are
  surfaced in `extract_components(what = "rank_diagnostics")` and
  `plot_components(plot_type = "ajive_diagnostic")` (forestgreen histogram
  and hline).
  **Argument order note:** `n_perm_samples` is positioned after `joint_rank`
  in the signature, preserving the original positional slot of `joint_rank`.

## New vignettes

- Added `vignettes/function_gallery.Rmd`: comprehensive API gallery covering
  all ~38 exported functions with simulated examples, plots, and
  interpretation notes.

---

# rajiveplus 0.1.0

## New features

- Added `rank_features()` with four modes: `top_loadings`, `contribution`,
  `overlap`, and `significant`.
- Added `summarize_components()` with unified summary modes:
  `components`, `variance`, `significance`, `associations`, `stability`, and
  `diagnostics`.
- Added remaining non-benchmark convenience modules:
  `export_results()`, `rajive_report()`, and compatibility wrappers for
  association and stability workflows.
- Extended `extract_components()` and `plot_components()` with additional
  unified modes for scores/loadings/variance/significance extraction and
  non-diagnostic plotting paths.
- Added compatibility wrappers: `get_top_loadings()`,
  `get_feature_contributions()`, `compare_feature_sets_across_blocks()`, and
  `summarize_significant_vars()`.
- Added tests covering return schema, ranking behavior, overlap summaries,
  significance summaries, and wrapper delegation.

- Added `autoplot.rajive()` and `autoplot.jackstraw_rajive()` S3 methods that
  delegate to `plot_components()` for seamless ggplot2 integration.
- Added `fortify.rajive()` and `fortify.jackstraw_rajive()` S3 methods that
  return tidy long-format data frames via `extract_components()`.
- Fixed `plot_components(plot_type = "pairs")` to degrade gracefully to a
  single-component scatter when only one joint component is available.

## Documentation

- Added vignette usage example for `rank_features()` in the CLL application
  workflow.
