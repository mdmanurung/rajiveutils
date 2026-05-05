# Statistical Audits

## Scope

This document records statistical audits of the current `rajiveplus` methods and
the planned interpretation / inference API.

Document boundary:
- This file owns methodological guardrails and interpretation constraints.
- `VisualizationPlan.md` owns API sequencing and function contracts.
- `PROGRESS.md` owns implementation/validation history once a guardrail is
  actioned in code.

## Cross-File Integration

- Loading-stability guardrails apply to `assess_stability(target = "loadings")`
  and `bootstrap_loading_stability()` in `VisualizationPlan.md`.
- Association-inference warnings apply to `associate_components()` and related
  score-association helpers.
- Survival split warnings apply to survival modes and any plotting/reporting
  wrapper that discretizes scores.
- Rank-selection wording constraints apply to diagnostic helpers and reporting
  text built around AJIVE thresholds.

## Critical Findings

### 1. Jackstraw Multiple-Testing Correction Must Be Global

Current implementation in `R/jackstraw.R` previously applied correction within
each block/component slice, and the Bonferroni threshold was computed as a
block-local quantity.

Statistical consequence:
- this is not a global family-wise correction across the full experiment
- block-level thresholds can overstate power relative to a global test family

Implemented action:
- default correction changed to global `BH`
- `BH` is now applied once across all block-component-feature tests
- `bonferroni` now divides by the total number of tests across all blocks and
  all joint components
- `jackstraw_rajive()` now stores `n_tests` as an attribute for downstream use

### 2. Current Jackstraw Null Is A Simplified Approximation

The current null generation permutes the score vector while holding the observed
decomposition fixed. This is a computational shortcut, not a full jackstraw
scheme that re-estimates latent structure after feature masking.

Statistical implication:
- results should be described as approximate jackstraw-style inference
- documentation should not imply full exact equivalence to the original
  jackstraw construction

Implemented action:
- `jackstraw_rajive()` documentation now states that the method uses a
  score-permutation shortcut
- downstream documentation should avoid language that implies full
  decomposition re-estimation

### 3. Procrustes Alignment Is Mandatory For Loading Stability

Any future bootstrap loading stability function is statistically invalid unless
bootstrap loadings are aligned to a reference solution up to sign and rotation.

Statistical implication:
- raw bootstrap loading variation is dominated by rotational indeterminacy
- unaligned variability summaries are not interpretable

Implemented action (2026-05-05 Session 31):
- implemented orthogonal Procrustes alignment in
  `assess_stability(target = "loadings")`
- loading-stability summaries (`mean_loading`, `sd_loading`,
  `cos_similarity`) are computed after alignment

### 4. Planned Association Helpers Need Explicit Inferential Warnings

Planned functions such as `associate_components()` will treat estimated AJIVE
scores as predictors or summaries. These scores carry estimation error.

Statistical implication:
- effect sizes are attenuated
- naive p-values can be misinterpreted if score estimation uncertainty is ignored

Implemented action (2026-05-05 Session 31):
- `associate_components()` now documents and emits an inferential warning that
  component-score association is post-decomposition and that score estimation
  error is not propagated into p-values

### 5. Planned Survival Split Helpers Need A Bias Warning

Any future survival helper using score splitting such as median/tertile groups
must state that split-based inference is data-adaptive and may be
anti-conservative.

Implemented action (2026-05-05 Session 31):
- `associate_components()` now warns that split-based survival inference is
  data-adaptive and may be anti-conservative
- documentation and runtime messaging prefer continuous-score Cox
  (`split = "none"`) as primary inference and labels split modes as descriptive

## Important But Non-Critical Findings

### 6. Rank Selection Thresholds Should Be Described More Precisely

The AJIVE rank selection combines:
- a Wedin-based threshold
- a random-direction threshold
- the rule `max(wedin, random_direction)`

This is a practical conservative rule, but it should not be described as having
formal FDR or FWER control for rank selection.

Implemented action (2026-05-05 Session 31):
- `plot_components(plot_type = "rank_threshold")` documentation now describes
  the threshold as a practical heuristic/conservative rule and explicitly
  avoids FWER/FDR overclaims

### 7. Variance Explained Language Should Be Tightened

For AJIVE, "variance explained" should distinguish:
- joint variance explained
- block-specific individual variance explained after joint signal removal

Implemented action (2026-05-05 Session 31):
- `plot_components()` documentation now sets language guardrails for future
  variance outputs: joint vs block-specific individual variance after joint
  removal must be reported distinctly

### 8. Jackstraw Significance Means Association With Estimated Scores

"Significant features" in jackstraw should be interpreted as features associated
with the estimated joint component scores, not necessarily causal drivers.

Implemented action (2026-05-05 Session 31):
- `jackstraw_rajive()` documentation now uses association-first wording and
  explicitly warns against causal overinterpretation

## Implementation Priority

### Completed (as of 2026-05-05 Session 27)
- global `BH` as default in `jackstraw_rajive()` — done; see Finding #1
- global Bonferroni (total test count) if `bonferroni` requested — done; see Finding #1
- documentation note on simplified score-permutation null — done; see Finding #2

### Completed in this session (2026-05-05 Session 31)
- **Finding #3**: Procrustes alignment added in `assess_stability(target = "loadings")`
- **Finding #4**: inferential warning implemented in `associate_components()`
- **Finding #5**: survival split bias warning implemented in `associate_components()`
- **Finding #6**: rank-threshold wording tightened as heuristic/non-FWER/non-FDR
- **Finding #7**: variance-language guardrail added for future variance plots
- **Finding #8**: jackstraw significance wording changed to association (non-causal)

### Follow-up still pending
- extend the same warning semantics to future compatibility wrappers/aliases
  once those APIs are added (for example, survival-specific convenience wrappers)