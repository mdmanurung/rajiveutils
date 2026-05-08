# rajiveplus (development)

## Bug fixes

- Fixed pkgdown site build failure: added `microbiome_application` vignette
  to the `articles` section of `_pkgdown.yml` so pkgdown no longer aborts
  with "vignette missing from index".
- Fixed pkgdown site build failure: marked all code chunks in
  `vignettes/benchmarking.Rmd` as `eval=FALSE` to prevent CI failures caused
  by the external `RaJIVE` GitHub install and missing cached `.rds` result
  files.  The benchmarking vignette now renders as a display-only reference;
  see the note at the top of the vignette for instructions on reproducing
  results locally.

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
