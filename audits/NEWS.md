# rajiveplus (development)

## Performance

- `jackstraw_rajive()` null F-statistic generation now reuses centered block
  rows and computes simple-regression F-statistics from sufficient statistics.
  The refactor preserves the previous RNG/result contract within `1e-10` in
  old-vs-new tests while reducing elapsed time and allocations in the
  PERF-002 benchmark.

## Bug fixes

- `extract_components(..., what = "variance", format = "long")` and
  `fortify.rajive(..., what = "variance")` now return tidy columns
  `block`, `component`, and `proportion`; the historical wide
  `Joint`/`Indiv`/`Resid` output remains the default for `format = "wide"`.

- Categorical and batch association helpers now accept only `NULL` or
  `"kruskal"` and report `"kruskal"` truthfully.  Unsupported labels such as
  `"anova"` now error with class `rajiveplus_invalid_input` instead of being
  reported without changing the underlying test.

- `assess_stability(method = "permutation")` now errors explicitly with class
  `rajiveplus_invalid_input`; permutation stability and `n_perm` remain
  reserved for a future implementation.

- `RobRSVD.all(nrank <= 0)` now returns `numeric(0)` singular values and
  zero-column `u`/`v` matrices before entering the C++ backend.

- `RobRSVD_all_cpp()` now zero-initialises the `d`/`u`/`v` accumulators
  and trims them to the number of components actually fitted before
  returning.  Previously, on the rare early-break path (when
  `arma::svd_econ` failed or returned an empty singular-value vector),
  trailing entries were uninitialised.  No realistic call exercised the
  faulty path, but the fix removes the latent UB.

- Native missing-data diagnostics and reconstructions now work with
  `full = FALSE`. Internal reconstruction falls back to `u/d/v` component
  fields when component `$full` matrices are intentionally omitted, while
  `get_block_matrix()` keeps its existing `NA` return for omitted full
  matrices.

- Native missing-data fits now store `fit$missing$estimability` and
  `fit$missing$reconstruction_provenance` during fit construction.
  `get_reconstructed_blocks()` attaches matching provenance as a
  `reconstruction_provenance` attribute on returned block lists.

- Native automatic rank diagnostics now score candidate joint ranks by
  holding out observed block/sample rows when possible, rather than by
  reusing the fitted observed-entry residual objective. This fixes the
  previous degeneracy where `joint_rank = NA` could always select the lowest
  candidate, including rank 0, because smaller joint ranks left more
  block-specific individual rank for in-sample reconstruction.

- Native missing-data fits now warn with class
  `rajiveplus_saturated_signal_rank` when any `initial_signal_ranks` entry is
  at or above that block's maximum matrix rank. Such settings can nearly
  interpolate observed entries and make residual variance look artificially
  close to zero.

- `rajive_missing_control()` gains `rank_repeats`, `svd_shrinkage`, and
  `svd_shrinkage_coeff`.
  Native automatic rank diagnostics now average held-out prediction error over
  repeated holdout splits by default. `svd_shrinkage` supports fixed
  soft-thresholding or `svd_shrinkage = "missmda"`, which shrinks retained
  singular values using a discarded-singular-value noise estimate inspired by
  missMDA's regularized imputation algorithms.

## Statistical changes

- `get_random_direction_bound_robustH()` now uses classical `base::svd()`
  to build the random-direction null instead of the robust M-estimator.
  The null is generated from i.i.d. Gaussian draws, where the M-estimator
  only adds Monte-Carlo noise without removing bias; the classical SVD
  matches the AJIVE reference implementation.  All other pipeline steps
  continue to use the robust SVD.  See
  `audits/2026-05-09-rajiveplus-vs-rajive-parity.md` for the rationale.

## New features

- Added native incomplete-data RaJIVE support behind
  `Rajive(..., missing = "native")`. The complete-data default remains
  `missing = "error"`. Native mode accepts an observed-entry `mask`, excludes
  missing cells from preprocessing and fitted objectives, stores
  `rajive_incomplete` metadata, and exposes `get_missingness_info()`,
  `get_estimability()`, `get_reconstructed_blocks()`,
  `get_missing_diagnostics()`, `plot_missingness()`,
  `diagnose_missing_ranks()`, and `get_missing_uncertainty()`.
  Reconstructions are explicitly derived outputs: observed values are never
  overwritten in the fit, and block-specific individual signal for an entirely
  missing sample-block row is labelled `not_identifiable`. Bootstrap
  uncertainty is a parametric residual bootstrap on observed entries;
  sensitivity reporting records requested assumptions but is not yet a full
  sensitivity analysis; and `diagnose_missing_ranks()` reports observed-entry
  prediction error and support diagnostics rather than null-model surrogates.

- Native missing-data fits now treat `joint_rank = NA` as automatic native
  rank selection, equivalent to `joint_rank = "native_cv"`. When no explicit
  `rank_candidates` are supplied, candidates default to
  `0:min(initial_signal_ranks)`, including the no-joint-signal rank. Fixed
  numeric `joint_rank` values remain fixed-rank fits and do not run automatic
  diagnostics.

- Added bootstrap inference helpers for the v0.2.0 roadmap:
  `rajive_ci()` computes percentile/basic/BCa intervals for joint loadings,
  variance explained, and joint rank, while rejecting BCa for sample-specific
  scores where the delete-one jackknife target is undefined.  The internal
  bootstrap engine supports observation, cluster, and stratified cluster
  resampling and can return aligned loading/score replicates for downstream
  reuse.

- Added `joint_variance_partition()` for feature-level joint/individual/
  residual sum-of-squares decomposition.  The function requires the original
  `blocks` and computes residuals explicitly as `X - J - I` instead of
  trusting the optional noise slot.

- `assess_stability()` can now attach bootstrap replicate arrays via
  `return_replicates = TRUE`, and accepts cluster/strata resampling arguments.
  Existing summaries are preserved by default.

- `associate_components()` gains `propagate_uncertainty = "bootstrap"`,
  adding stability, effect-size quantile, and median-p-value summaries while
  preserving the original point-estimate `p_value` output.

- `extract_components(..., ci = ...)`, `plot_components(plot_type =
  "loading_ci")`, and `autoplot(..., plot_type = "loading_ci")` provide a
  user-facing path from `rajive_ci()` output to loading-interval tables and
  plots.

- `jackstraw_rajive()` gains three new arguments for posterior inclusion
  probabilities: `pip` (default `FALSE`), `pip_pi0`, and `pip_group`
  (`"component"`, `"block_component"`, or `"pooled"`).  When `pip = TRUE`,
  each block/component result list includes a `$pip` vector of per-feature
  PIPs computed as `1 - qvalue::lfdr(p_values)`.  Requires the
  **qvalue** Bioconductor package (Suggests).
  Reference: Chung NC (2020). *Bioinformatics*, 36(10):3107–3114.

- `Rajive()` gains a `seed` argument (default `NA`).  When non-`NA`,
  `set.seed(seed)` is called before any random sampling, making results
  fully reproducible without setting the seed externally.

- `jackstraw_rajive()` gains a `pool` argument (default `"block"`).  When
  `"block"`, null F-statistics from all joint components are pooled within
  each block before computing empirical p-values, giving a
  $d_k \times (J \cdot n_{null})$ null pool per block.  Set `pool = "global"`
  to restore the original per-component behaviour.

- `Rajive()` now computes both the permutation bound and the random-direction
  bound independently when both `n_perm_samples` and `n_rand_dir_samples` are
  supplied (previously the permutation bound silently replaced the
  random-direction bound).  The final threshold is still `max(wedin,
  rand_dir, perm)`.

- `assess_stability(target = "components")` now attaches a `rank_match_rate`
  attribute to the returned data frame: the fraction of successful bootstrap
  iterations in which the bootstrap joint rank equalled the reference rank.
  Values below 1 indicate rank instability across resamples.

## Changes in defaults

- `jackstraw_rajive()`: default `correction` reverted from `"BY"` back to
  `"BH"` (Benjamini--Hochberg).  BH is the standard FDR procedure and
  sufficient for most use cases; switch to `correction = "BY"` explicitly
  when neighbouring features are strongly correlated.

- `jackstraw_rajive()`: default `pip_group` changed from `"component"` to
  `"pooled"`.  A single pooled `lfdr` call over all tests gives a more
  stable `π₀` estimate when individual blocks are small.  Use
  `pip_group = "component"` to restore the previous behaviour.

## Bug fixes / behaviour changes

- Degenerate columns (near-zero variance) in a `Rajive()` input block are now
  **dropped automatically** with a warning (`class =
  "rajiveplus_degenerate_block"`) instead of aborting with an error.  This
  allows analysis to continue when a small fraction of constant features are
  present (e.g. batch-invariant housekeeping genes).

- `truncate_svd(rank = 0)` now returns zero-column `u` and `v` matrices
  instead of `n x 1` zero placeholders.  This keeps internal rank bookkeeping
  consistent for zero-rank decompositions.

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

- Added `vignettes/inference.Rmd`: simulated-data inference walkthrough that
  exercises feature-level variance partition, loading CIs, BCa rank intervals,
  loading-CI plotting, and score-uncertainty propagation.  A fast render check
  is included in the test suite when the source vignette is available.

- Added `vignettes/function_gallery.Rmd`: comprehensive API gallery covering
  all ~38 exported functions with simulated examples, plots, and
  interpretation notes.

- Added `vignettes/native_missing_union.Rmd`: a portable synthetic
  BMV-like example showing union-sample alignment across partially
  overlapping blocks, native missing-data fitting with `full = FALSE`,
  missingness diagnostics, estimability labels, and reconstructed block
  provenance.

- Added `vignettes/microbiome_application.Rmd`: end-to-end RaJIVE workflow
  on a multi-kingdom gut microbiome dataset (Haak et al. 2021; bacteria,
  fungi, viruses). Demonstrates rank selection, variance decomposition,
  joint/individual score inspection, association testing, and bootstrap
  stability.

## Other changes

- The heavy benchmarking artifact moved from `vignettes/benchmarking_heavy.Rmd`
  to `inst/benchmarks/benchmarking_heavy.Rmd`.  It remains SLURM-only and
  writes/reads caches explicitly under repo `vignettes/data`.

- Added local artifact `inst/analyses/bmv_native_missing_union.Rmd` plus
  `logs/slurm_render_bmv_native_missing.sh` to render the BMV union-sample
  native missing-data analysis outside package checks. The artifact restores
  `NA` support from the BMV preprocessing cache and keeps full BMV fitting
  behind `RUN_FULL_BMV_NATIVE=1`.

- Added audit regression tests for variance fortification, association method
  validation, stability method validation, zero-rank robust SVD, and
  pkgdown/vignette metadata.

- `data_heatmap()` is now exported. The internal `geom_raster()` mapping now
  uses `scale_y_discrete()` / `scale_x_discrete()` (previously
  `scale_y_continuous()` / `scale_x_continuous()`), which correctly matches
  the factor-encoded axes and eliminates ggplot2 warnings.

- `get_perm_bound_robustH()` is now marked `@keywords internal`; it no longer
  appears in the pkgdown reference index.

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
