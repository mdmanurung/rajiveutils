# Scientific & Statistical Audit — `rajiveplus` (full package)

**Date:** 2026-05-08
**Auditor mode:** read-only
**Scope:** primary analysis paths end-to-end (`Rajive()` → joint-rank selection → final decomposition; `jackstraw_rajive()`; `associate_components()`; `assess_stability()`; supporting numerics in `RobRSVD_*` and simulation helpers)
**Files inspected:** [R/Rajive.R](../R/Rajive.R), [R/Rajive_helpfunctions.R](../R/Rajive_helpfunctions.R), [R/jackstraw.R](../R/jackstraw.R), [R/visualization.R](../R/visualization.R), [R/RobustSVD.R](../R/RobustSVD.R), [R/simulation.functions.R](../R/simulation.functions.R), [src/RobustSVD.cpp](../src/RobustSVD.cpp), [tests/testthat/test-jackstraw.R](../tests/testthat/test-jackstraw.R), [tests/testthat/test-rajive-perm.R](../tests/testthat/test-rajive-perm.R)

> **Read-only invariant honored.** No source, test, doc, or config file was modified. Only this report was written. No `devtools::document()` / install / commit was run.

---

## 1. Executive Summary

- **Overall verdict — partially defensible, with several uncalibrated claims.** The core RaJIVE machinery (block SVD → stacked signal-score SVD → identifiability filter → final J/I/E decomposition) is implemented coherently and matches the published RaJIVE algorithm. However, several quantitative claims (joint-rank thresholds, jackstraw p-values, bootstrap "stability" of components) rest on heuristic null distributions or shortcut estimators that are **not formally calibrated** under any explicit error-rate guarantee and are only partly disclosed to the user.
- **High-risk: `sim_dist()` silently fails for `dist.type = 3` (exponential).** The dispatch has no `num == 3` branch; the variable `dist` is left undefined and the function errors with an opaque message — but the documentation and the public `ajive.data.sim()` parameter advertise this as a supported option. Any simulation-based test of robustness against heavy-tailed noise is currently impossible via the public API. ([R/simulation.functions.R:80-99](../R/simulation.functions.R#L80-L99))
- **High-risk: `get_sv_threshold()` returns `NA` when `rank == length(singular_values)`,** because it indexes `sv[rank+1]`. The downstream comparison `d > NA → NA`, then `sum(NA, ...) = NA`, and `truncate_svd(rank = NA)` will error. When `initial_signal_ranks[k]` equals the full SVD rank — a normal case for small / square blocks — the pipeline silently degrades or aborts. ([R/Rajive_helpfunctions.R:25-29](../R/Rajive_helpfunctions.R#L25-L29), [R/Rajive.R:115-118](../R/Rajive.R#L115-L118))
- **High-risk: parallel RNG is not configured.** `mclapply()` and `foreach::%dopar%` paths in `Rajive()`, `get_wedin_bound_samples()`, `get_random_direction_bound_robustH()`, `get_perm_bound_robustH()`, and `assess_stability()` all draw random numbers with `num_cores > 1`, but no `RNGkind("L'Ecuyer-CMRG")` / `mc.set.seed = TRUE` / `doRNG::%dorng%` is used. **Any user that sets `set.seed()` and then runs with `num_cores > 1` will get non-reproducible results.** ([R/Rajive_helpfunctions.R:106-137, 150-176, 196-220](../R/Rajive_helpfunctions.R#L106-L220), [R/Rajive.R:62-86](../R/Rajive.R#L62-L86))
- **Medium-risk: jackstraw is a permute-the-scores shortcut, not a proper jackstraw.** The implementation permutes the *score vector* against fixed feature rows rather than refitting the latent decomposition after permuting `s` synthetic features (Chung & Storey 2015). The disclaimer in the docstring is appropriate, but the test suite never compares calibration against a real jackstraw or against the analytic F null. ([R/jackstraw.R:135-166, 244-271](../R/jackstraw.R#L135-L271))
- **Medium-risk: joint-rank threshold combination has no formal calibration.** `overall_sv_sq_threshold = max(wedin_5%, rand_dir_95%, perm_95%)` (with `na.rm = TRUE`) is a heuristic; the per-test FWER/FDR after this max-rule is not bounded. The Wedin-bound resampling itself uses bootstrap-with-replacement of perp-basis columns rather than uniform sampling on a sphere. ([R/Rajive.R:201-226](../R/Rajive.R#L201-L226), [R/Rajive_helpfunctions.R:127-138](../R/Rajive_helpfunctions.R#L127-L138))
- **Medium-risk: `assess_stability(target = "components")` never aligns bootstrap components to the reference and never takes absolute correlation.** Sign flips and component re-ordering will systematically depress reported "mean_correlation" and may incorrectly suggest instability. Procrustes alignment is applied to the loadings branch but not the components branch. ([R/visualization.R:1535-1570](../R/visualization.R#L1535-L1570))
- **Medium-risk: identifiability filter index recomputation looks suspicious.** The drop loop appends `j` and breaks the inner-`k` loop, but the loop variable `j` always refers to the column index of the *full* `joint_scores` matrix from `M_svd[['u']]`. This is correct, but the post-loop `joint_rank <- length(to_keep)` does not log which columns were originally dropped to the user via any user-facing summary; only `rank_sel_results[['identif_dropped']]` is set. Print/summary methods do not surface it.

Sections 2–5 below give the verifiable evidence and a runnable Verification Plan.

---

## 2. High-Risk Scientific / Statistical Issues (Definite Errors)

### 2.1 `sim_dist()` silently broken for `dist.type = 3`

- **Finding** — `dist.type = 3` ("exponential") is documented and exposed via `ajive.data.sim()` but `sim_dist()` has no `num == 3` branch; `dist` is undefined → `Error: object 'dist' not found`. Not caught by any test.
- **Location** — [R/simulation.functions.R:88-99](../R/simulation.functions.R#L88-L99); doc claim at [R/simulation.functions.R:9](../R/simulation.functions.R#L9) and [R/simulation.functions.R:74](../R/simulation.functions.R#L74).
- **Evidence** —
  ```r
  if (num == 1) {
    dist <- rnorm(n*p)
  }else if (num == 2) {
    dist <- runif(min=0, max=1, n=n*p)
  }
  out <- matrix(dist, nrow = n, ncol = p)
  ```
  No `else if (num == 3)` block; no `else stop(...)` either.
- **Consequence** — Any user attempting simulation under heavy-tailed noise (the documented use case) gets a confusing R error rather than the advertised exponential block. Robustness claims for RaJIVE that depend on the simulation harness cannot be validated using the package's own `ajive.data.sim()`.

### 2.2 `get_sv_threshold()` returns `NA` at the boundary `rank == length(d)`

- **Finding** — Threshold = `0.5 * (d[rank] + d[rank+1])`; when `rank == length(d)`, `d[rank+1]` is `NA`, threshold becomes `NA`. This is then passed to `sum(d > sv_threshold)` (NA propagation) and `truncate_svd(rank = NA)`.
- **Location** — [R/Rajive_helpfunctions.R:25-29](../R/Rajive_helpfunctions.R#L25-L29) (and a second copy in [R/Rajive.R:114-118](../R/Rajive.R#L114-L118) — note the duplication).
- **Evidence** — The fast path in `Rajive()` deliberately requests only `signal_rank + 1` singular values:
  ```r
  svd_ranks <- pmin(initial_signal_ranks + 1L, max_rank_per_block)
  block_svd <- parallel::mclapply(seq_along(blocks), function(k)
    get_svd_robustH(blocks[[k]], rank = svd_ranks[k]), mc.cores = num_cores)
  ```
  When `initial_signal_ranks[k] >= min(dim(X_k)) - 1`, `svd_ranks[k] == min(dim(X_k))`, so `singular_values` has exactly `rank` entries and `d[rank+1] = NA`.
- **Consequence** — Silent NA propagation into `sv_thresholds`, then either:
  (a) the inner `if(sv < sv_thresholds[[k]])` (Rajive.R line ~280) compares to NA → `NA < x` → `NA` → `if(NA)` errors with *"missing value where TRUE/FALSE needed"*; or
  (b) downstream individual rank selection (`indiv_rank <- sum(d > sv_threshold)`) returns `NA`, and `truncate_svd` indexes `1:NA` → error.
  Either way the user sees a spurious failure instead of a clear "rank too close to dim(X)" message.

### 2.3 Parallel RNG: results are non-reproducible with `num_cores > 1`

- **Finding** — Wedin, random-direction, and permutation bound samplers all parallelize via `mclapply()` / `foreach::%dopar%` over `doParallel`. Neither `RNGkind("L'Ecuyer-CMRG")` nor `mc.set.seed = TRUE` is used; `foreach` is not switched to `doRNG::%dorng%`. With `num_cores > 1`, the random samples drawn inside each worker are not seeded reproducibly even when the user calls `set.seed()` before `Rajive()`.
- **Locations** —
  - `mclapply` blocks: [R/Rajive.R:67, 76, 82](../R/Rajive.R#L67-L82).
  - `%dopar%` blocks: [R/Rajive_helpfunctions.R:128](../R/Rajive_helpfunctions.R#L128), [R/Rajive_helpfunctions.R:165](../R/Rajive_helpfunctions.R#L165), [R/Rajive_helpfunctions.R:215](../R/Rajive_helpfunctions.R#L215).
  - `assess_stability()` parallel-safe wrapper is single-threaded and uses `tryCatch()` around the bootstrap fit, so reproducibility there is okay *only* when `num_cores = 1` is the default *and* `Rajive()` itself is not threaded — but `assess_stability()` still calls `Rajive(b_list, initial_signal_ranks)` with the default `num_cores = 1L`, so this particular call is reproducible. The exposure is on direct user calls to `Rajive(..., num_cores = K)`.
- **Evidence** — `foreach::foreach(...) %dopar% { ... rnorm(n_obs * l, ...) ... }` (random-direction bound) and `... sample.int(n_obs) ...` (perm bound) are dispatched to worker processes that inherit the parent's RNG state once but then advance independently and irreproducibly under fork-based parallelism without an L'Ecuyer stream.
- **Consequence** — Any peer reviewer or downstream user who sets a seed and tries to re-run the analysis will get *different* joint-rank thresholds and possibly different joint ranks. This silently breaks reproducibility — a load-bearing scientific property for any inferential method.

---

## 3. Medium-Risk Issues (Likely Problems)

### 3.1 Jackstraw is a score-permutation shortcut, not a refit-after-mask jackstraw

- **Finding** — The Chung & Storey (2015) jackstraw permutes a small number of randomly-chosen feature rows and **refits the latent factor model** so the resulting null F-statistics absorb the impact of permuted rows on the estimated scores. This implementation instead fixes `joint_scores` (estimated once on the full data) and permutes the score vector for null draws — i.e. it tests `H0: cor(feature, joint_score_j) = 0` with the joint scores treated as an exogenous covariate, ignoring the data-dependence of the scores on the feature itself.
- **Location** — [R/jackstraw.R:135-166](../R/jackstraw.R#L135-L166) (`generate_null_f_stats`), [R/jackstraw.R:244-271](../R/jackstraw.R#L244-L271) (main loop). The docstring at [R/jackstraw.R:170-186](../R/jackstraw.R#L170-L186) discloses this as "approximate jackstraw-style".
- **Evidence** — `null_f[, s] <- ols_f_stat_matrix(rows, scores_perm)` uses `scores_perm <- sample(joint_comp_scores)` rather than refitting Rajive on a partially-permuted block matrix.
- **Consequence** — When a feature has a strong loading on the true joint signal, the joint score estimated *with that feature included* is partially driven by the feature itself. Permuting the score breaks the marginal feature↔score association but does **not** undo the contribution of the feature to the score. Net effect: **anti-conservative p-values**, biased toward declaring more features significant than the proper jackstraw would. Magnitude is unknown; depends on `d_k / n` and feature-to-signal SNR.
- **What would resolve the uncertainty** — Calibration simulation under the global null (`H0`: random Gaussian blocks with zero true joint structure). The empirical CDF of `p_values` should be `~ U(0,1)`. Run the simulation in Section 5.1.

### 3.2 Empirical p-value uses strict `>` rather than `≥`

- **Finding** — `compute_empirical_pvalues()` uses `n_ge = N - findInterval(f_obs, sorted_pool)` which equals `#{pool > f_obs}`, but the Phipson–Smyth formula counts `#{pool ≥ f_obs}`. The doc claims `≥`.
- **Location** — [R/jackstraw.R:97-110](../R/jackstraw.R#L97-L110).
- **Evidence** — `findInterval(x, vec)` (default `left.open = FALSE`) returns `max{i : vec[i] ≤ x}`, so `N - n_le = #{i : vec[i] > x}`, not `≥`.
- **Consequence** — For continuous F-statistics with no exact ties this is negligible. With discrete or rounded null pools (e.g., constant-feature degeneracies, very small `n`), p-values are slightly too small (anti-conservative) by 1/(N+1) per tie.
- **What would resolve the uncertainty** — Add `n_ge = N - findInterval(f_obs, sorted_pool, left.open = TRUE)` and unit-test against `mean(pool >= f_obs)`.

### 3.3 Joint-rank threshold = `max(wedin_5%, rand_dir_95%, perm_95%)` has no calibrated error rate

- **Finding** — Three thresholds with different statistical interpretations are combined by `max(..., na.rm = TRUE)`, then applied component-by-component to the squared singular values. The docstring at [R/visualization.R:196-200](../R/visualization.R#L196-L200) acknowledges "does not carry formal FWER or FDR" but the user-facing `Rajive()` documentation does not.
- **Location** — [R/Rajive.R:201-226](../R/Rajive.R#L201-L226).
- **Evidence** — The Wedin samples target the 5th-percentile *lower* bound on the leading squared SV under the alternative; the random-direction and permutation samples target the 95th-percentile *upper* bound on the leading squared SV under the null. Taking `max` produces a conservative upper threshold for the *first* component, but the same threshold is applied to *every* component of `M_svd[['d']]` via `sum(d^2 > t)`. Multiple-testing across components is not corrected.
- **Consequence** — Bias toward too-large joint ranks when the true rank is `> 1`, because the null bounds are calibrated for the leading SV only. Magnitude is data-dependent.
- **What would resolve the uncertainty** — Calibration simulation under `rankJ = 0` (Section 5.2): the empirical false-positive rate of "joint_rank > 0" should match the nominal 5% if the null bound is calibrated; observed inflation as `K` and per-block dimensions grow would confirm this concern.

### 3.4 Wedin bound uses bootstrap of perp-basis columns, not uniform sampling on a sphere

- **Finding** — The classical Wedin sin-θ resampling samples Haar-uniform directions from the perp subspace; this implementation does `sampled_col_index <- sample.int(ncol(perp_basis), size = ncol(perp_basis), replace = TRUE)` and uses the bootstrapped basis matrix directly. The resulting matrix is *not* an orthonormal basis (columns are duplicated and non-orthogonal).
- **Location** — [R/Rajive_helpfunctions.R:127-138](../R/Rajive_helpfunctions.R#L127-L138).
- **Evidence** — `perp_resampled <- perp_basis[ , sampled_col_index]` then `norm(X %*% perp_resampled, type = '2')` — the operator-2 norm of `X * Q_boot` is *not* an unbiased estimate of `||X * Q_haar||_2` because `Q_boot` has Frobenius norm `sqrt(rank)` but is not orthonormal; the computed norm is biased upward by the column-overlap structure.
- **Consequence** — Wedin samples are slightly *larger* than the proper Haar-Wedin samples → 5th-percentile threshold is *higher* → since this enters as `max(wedin, rand_dir, perm)`, the joint rank is biased *downward* (too conservative) when Wedin dominates. This may explain why the rand-dir threshold typically dominates in practice.
- **What would resolve the uncertainty** — Replace `perp_basis[, sample.int(...)]` with a Haar-uniform draw (`qr.Q(qr(matrix(rnorm(p*r), p, r)))` projected through `perp_basis`) and compare distributions on a controlled simulation.

### 3.5 `assess_stability(target = "components")` ignores sign / order ambiguity of singular vectors

- **Finding** — Bootstrap correlation between `ref[idx, j]` and `bs[, j]` is computed as a signed correlation with the j-th bootstrap component matched to the j-th reference component by index. Singular vectors have arbitrary sign; equal singular values have arbitrary rotation; bootstrap re-fits routinely permute and flip signs.
- **Location** — [R/visualization.R:1551-1570](../R/visualization.R#L1551-L1570).
- **Evidence** —
  ```r
  for (j in seq_len(n_use)) {
    cor_mat[b, j] <- suppressWarnings(stats::cor(ref[idx, j], bs[, j],
                                                  use = "pairwise.complete.obs"))
  }
  ```
  No `abs()`, no Procrustes/Hungarian alignment. Compare the loadings branch ([R/visualization.R:1505-1510](../R/visualization.R#L1505-L1510)), which *does* call `.procrustes_align`.
- **Consequence** — `mean_correlation` is systematically biased toward 0 when sign flips occur (exactly the regime where the user expects "stable" components). Users will see a stable component reported as `mean_correlation ≈ 0` and incorrectly conclude the decomposition is unstable.
- **What would resolve the uncertainty** — Either `cor_mat[b, j] <- abs(cor(...))` or align bootstrap scores to reference via Procrustes / matching rotation before correlation.

### 3.6 `associate_components(mode = "survival")` coerces status with `as.logical(s_vec)`

- **Finding** — `survival::Surv(as.numeric(t_vec), as.logical(s_vec))` works for 0/1 numeric coding; for character `"alive"/"dead"` it silently coerces to `NA`, producing an all-censored Surv object. No validation upstream.
- **Location** — [R/visualization.R:1320-1330](../R/visualization.R#L1320-L1330).
- **Evidence** — `as.logical(c("alive","dead"))` → `c(NA, NA)` in R; no warning is emitted by `Surv()` either when many statuses are NA (it sets event = NA).
- **Consequence** — Cox model is fit on an all-censored or mostly-censored dataset; coefficient and p-value are meaningless but no error is raised.
- **What would resolve the uncertainty** — Add an input check (`stop()` if `is.character(s_vec)`) — but this is a fix, not part of this read-only audit; flagged for later.

### 3.7 `RobRSVD1_cpp` uses `arma::solve(uterm1, …)` without rank or conditioning safeguards

- **Finding** — Inside the M-estimator iteration the linear systems `solve(uterm1) %*% uterm2` and `solve(vterm1) %*% vterm2` are formed from sums of weighted outer products. When `mysigma → 0` (median absolute residual is zero, e.g., perfectly low-rank input or all-zero block), `Wmat = huberk / |R/sigma|` blows up and `uterm1`/`vterm1` become rank-deficient or non-finite.
- **Location** — [src/RobustSVD.cpp:64-101](../src/RobustSVD.cpp#L64-L101).
- **Evidence** — The R-side guard at [R/Rajive.R:330-345](../R/Rajive.R#L330-L345) (`get_joint_decomposition_robustH` zero-rank path) was added precisely because the C++ M-estimator "fails to converge: `solve(): solution not found`" on a zero matrix. Other near-degenerate inputs (constant rows/columns, perfectly low-rank blocks) are still vulnerable.
- **Consequence** — On real datasets with feature blocks that contain near-zero variance columns or perfectly co-linear features, `RobRSVD_all_cpp` may throw `arma::solve` errors with no informative context. Users see a low-level Armadillo error instead of a domain-meaningful message.

### 3.8 Sample-size-adequacy / identifiability for joint rank vs. `n`

- **Finding** — `Rajive()` does not check that `joint_rank_estimate <= min(initial_signal_ranks)` or that `n` is large enough relative to `sum(initial_signal_ranks)` for the stacked signal-score matrix `M` (which is `n × sum(initial_signal_ranks)`) to admit a meaningful SVD. When `n < sum(initial_signal_ranks)`, `M` has rank at most `n` and many of its singular values are 0; the joint rank estimate is then determined entirely by which of these zeros squeak above the (possibly tiny) thresholds.
- **Location** — [R/Rajive.R:172-176](../R/Rajive.R#L172-L176).
- **Consequence** — In high-dimensional / small-`n` settings (typical in single-cell or microbiome — explicit target use cases per `PLANS.md`), the joint rank can be poorly identified without any warning.

---

## 4. Unclear Assumptions Requiring Human Domain Review

Phrased for a domain expert reviewer:

1. **Robustness target.** RaJIVE's robust SVD uses Huber M-estimation with `huberk = 1.345`. Is that calibration appropriate for the scale of *centered, log-transformed microbiome counts* (the documented use case in `microbiome_application.Rmd`), or should `huberk` be tuned per-block?
2. **Block normalization.** The package does not enforce per-block scaling. Is the practitioner expected to standardize blocks (Frobenius-norm rescaling à la JIVE) before calling `Rajive()`? The vignettes do center but do not Frobenius-rescale; this means a high-variance block dominates the stacked-signal SVD `M`. Document expectation explicitly.
3. **Permutation bound exchangeability.** `get_perm_bound_robustH()` permutes rows of each block's signal scores `U_k` independently. This null preserves marginal block structure but assumes samples are *exchangeable across permutations*. In paired/longitudinal designs (matched pre/post), this is the wrong null. Should the API accept a `strata` argument?
4. **What does "associated with the joint score" mean to your end user?** The jackstraw test, as implemented, regresses each feature on the *estimated* joint score (post-decomposition). A significant feature is one whose values correlate with the estimated score, not necessarily one that drives the latent variable. Is this distinction made clear in your manuscript / vignette?
5. **Survival-analysis use case.** `associate_components(mode = "survival", split = "median")` is data-adaptive (cutpoint at observed median). The current warning mentions anti-conservativeness but does not point to any correction (e.g., `maxstat`-style permutation p-value). Is the intended use exploratory only?

---

## 5. Verification Plan

For each item in §2–§4, a concrete, runnable check.

### 5.1 Calibration of jackstraw p-values under global null (medium #3.1, #3.2)

```r
# Pseudocode — DO NOT modify package; run from a scratch script.
library(rajiveplus)
set.seed(2026)
B  <- 200
ps <- replicate(B, {
  K  <- 2; n <- 60; pks <- c(40, 30)
  blk <- lapply(pks, function(p) matrix(rnorm(n*p), n, p))   # null: no joint structure
  fit <- Rajive(blk, initial_signal_ranks = c(5, 4),
                joint_rank = 1, num_cores = 1L)              # force joint_rank=1
  js  <- jackstraw_rajive(fit, blk, n_null = 50L, alpha = 0.05,
                          correction = "none")
  unlist(lapply(js, function(b) lapply(b, `[[`, "p_values")))
})
ks.test(ps, "punif")             # expected p > 0.1 under H0 if calibrated
mean(ps <= 0.05)                 # expected ~ 0.05 if calibrated
```
Failure mode confirms anti-conservativeness from §3.1.

### 5.2 Calibration of `joint_rank` selection under H0 (medium #3.3, #3.4)

```r
set.seed(2026)
B <- 200
ranks_under_H0 <- replicate(B, {
  blk <- lapply(c(50, 40, 30), function(p) matrix(rnorm(60*p), 60, p))
  fit <- Rajive(blk, initial_signal_ranks = c(6, 5, 4),
                n_wedin_samples = 200, n_rand_dir_samples = 200, num_cores = 1L)
  fit$joint_rank
})
mean(ranks_under_H0 > 0)         # nominal: 0.05; >> 0.05 confirms inflation
table(ranks_under_H0)
# Repeat with n_perm_samples = 200 in place of n_rand_dir_samples for §3.3 perm path.
```

### 5.3 Recovery under known signal (sanity)

```r
set.seed(2026)
sims <- replicate(50, {
  Y <- ajive.data.sim(K = 3, rankJ = 2, rankA = c(5, 4, 3),
                      n = 80, pks = c(60, 50, 40), dist.type = 1)
  fit <- Rajive(Y$sim_data, initial_signal_ranks = c(5, 4, 3))
  fit$joint_rank
})
mean(sims == 2)                  # expected close to 1
```

### 5.4 `sim_dist(num = 3)` reachability (high #2.1)

```r
testthat::expect_error(rajiveplus:::sim_dist(num = 3, n = 5, p = 5),
                       regexp = "exponential|object 'dist' not found")
# This test will pass today. Once num=3 is implemented, replace with positive expectation.
```

### 5.5 `get_sv_threshold()` boundary (high #2.2)

```r
testthat::expect_no_error({
  set.seed(1)
  X <- matrix(rnorm(20), 5, 4)                     # min(dim)=4
  fit <- rajiveplus:::get_svd_robustH(X, rank = 4) # rank == min(dim)
  rajiveplus:::get_sv_threshold(fit$d, rank = 4)   # currently returns NA
})
# Currently FAILS by producing NA, propagated downstream.
```

### 5.6 Parallel RNG reproducibility (high #2.3)

```r
Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = 40, pks = c(30, 20))
set.seed(7); a <- Rajive(Y$sim_data, c(5,4), n_rand_dir_samples = 100, num_cores = 2)
set.seed(7); b <- Rajive(Y$sim_data, c(5,4), n_rand_dir_samples = 100, num_cores = 2)
testthat::expect_equal(a$joint_rank_sel$rand_dir$rand_dir_samples,
                       b$joint_rank_sel$rand_dir$rand_dir_samples)
# Expected: FAILS today (samples differ). Will PASS once L'Ecuyer streams are wired in.
```

### 5.7 `assess_stability(components)` sign/order invariance (medium #3.5)

```r
# Property: bootstrapping with B=2 on identical data should yield
# |mean_correlation| close to 1 for true components.
set.seed(7)
Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = 60, pks = c(40, 30))
fit <- Rajive(Y$sim_data, c(5, 4))
stab <- assess_stability(fit, Y$sim_data, c(5, 4),
                         target = "components", B = 20L, sample_frac = 0.9)
# Today: many entries near 0 due to sign flips.
# Property check after fix: all(abs(stab$mean_correlation) > 0.7)
```

### 5.8 Wedin Haar vs bootstrap discrepancy (medium #3.4)

Compare empirical distribution of `wedin_samples` (current implementation) against a Haar-uniform reference computed by drawing `Q ~ qr.Q(qr(matrix(rnorm(p*r), p, r)))` and projecting through `perp_basis`. KS-test on the two samples; large `D` confirms bias.

### 5.9 Survival status coercion (medium #3.6)

```r
meta <- data.frame(t = runif(40, 1, 10), status = sample(c("alive","dead"), 40, TRUE))
fit  <- Rajive(ajive.data.sim(K=2, rankJ=2, rankA=c(5,4), n=40, pks=c(30,20))$sim_data,
               c(5,4))
res  <- associate_components(fit, meta, mode = "survival",
                             time_col = "t", status_col = "status",
                             variable = "t")
# Today: silently fits an all-censored Cox model. After validation: should stop() with
# a clear "status_col must be 0/1 numeric or logical" message.
```

### 5.10 Block normalization sensitivity (unclear #4.2)

Rescale one block by `1e3`; re-run Rajive with otherwise identical input; compare joint rank and joint_scores. If results are not equivariant under per-block scaling, document the expectation.

---

## Verification Standard / Trace Blocks

### Trace: `Rajive()` joint-rank selection path

```
user → Rajive(blocks, initial_signal_ranks, n_wedin_samples, n_rand_dir_samples,
              joint_rank, n_perm_samples, num_cores)            R/Rajive.R:50
   └─ get_svd_robustH(X_k)  (per block, mclapply)                R/Rajive.R:67/82
        └─ RobRSVD.all → RobRSVD_all_cpp                          R/RobustSVD.R:14
   └─ get_sv_threshold(d_k, rank=initial_signal_ranks[k])         R/Rajive.R:90
        ←  ⚠ NA when rank == length(d) (Finding 2.2)
   └─ get_joint_scores_robustH(...)                               R/Rajive.R:97
        ├─ stack signal scores → M; M_svd                          R/Rajive.R:170-176
        ├─ if !is.na(n_wedin_samples): get_wedin_bound_samples     R/Rajive.R:183-198
        │     └─ wedin_bound_resampling   (⚠ Finding 3.4)         R/Rajive_helpfunctions.R:120-138
        ├─ if !is.na(n_perm_samples): get_perm_bound_robustH       R/Rajive.R:206-216
        │     ⚠ parallel RNG not seeded (Finding 2.3)             R/Rajive_helpfunctions.R:196-220
        ├─ elif !is.na(n_rand_dir_samples): get_random_direction…  R/Rajive.R:218-227
        │     ⚠ parallel RNG not seeded; uses RobRSVD on N(0,I)
        ├─ overall_t = max(wedin_5%, rand_95%, perm_95%, na.rm=T)  R/Rajive.R:229-231
        │     ⚠ no calibrated FWER/FDR (Finding 3.3)
        └─ identifiability filter (L2 norm vs sv_thresholds)        R/Rajive.R:265-285
   └─ get_final_decomposition_robustH (per block, mapply)          R/Rajive.R:111
```

### Trace: `jackstraw_rajive()`

```
user → jackstraw_rajive(ajive_output, blocks, alpha, n_null, correction)  R/jackstraw.R:235
   └─ for k in K:  X_t = t(blocks[[k]])                                    R/jackstraw.R:259
        for j in joint_rank:
          f_obs  ← ols_f_stat_matrix(X_t, scores_j)                         R/jackstraw.R:266
          f_null ← generate_null_f_stats(X_t, scores_j, n_null)             R/jackstraw.R:267
                  (⚠ Finding 3.1: permute scores, not refit Rajive)
          p_vals ← compute_empirical_pvalues(f_obs, f_null)                 R/jackstraw.R:268
                  (⚠ Finding 3.2: '>' instead of '>=')
   └─ all_p_adj ← p.adjust(unlisted, method = correction)                   R/jackstraw.R:294-300
```

### Claim provenance

- "Tests pass" claims in `PROGRESS.md` (207 tests) verify *execution* and *consistency with previous behaviour*. They do **not** verify statistical calibration (no null-simulation tests, no recovery-under-known-signal property tests). Specifically:
  - [tests/testthat/test-jackstraw.R](../tests/testthat/test-jackstraw.R) covers `ols_f_stat_matrix` correctness (good — matched against `lm()`), but `compute_empirical_pvalues` tests use `runif`-generated null pools and check sign of the p-value, not uniformity under H0.
  - [tests/testthat/test-rajive-perm.R](../tests/testthat/test-rajive-perm.R) checks that the perm path *runs* and *replaces* rand_dir; it does not check that `joint_rank` under simulated H0 has the nominal Type-I rate.
- The audit-fix history in `PROGRESS.md` (sessions 31, 38, 39) addresses earlier audit findings #3–#8 mostly via documentation/warning messages, *not* via algorithmic correction. The Finding-#3 fix (L2 norm replacing operator-1 norm in the identifiability filter) is real and correctly applied.

---

## Pointer / Executive Summary for chat

Audit report written to [audits/2026-05-08-full-package.md](audits/2026-05-08-full-package.md).

Headline: core RaJIVE algebra is implemented coherently, but **(i)** `sim_dist(dist.type = 3)` is silently broken, **(ii)** `get_sv_threshold()` returns NA at the boundary `rank == length(d)`, **(iii)** parallel paths have no L'Ecuyer RNG so `set.seed()` does not reproduce results when `num_cores > 1`, **(iv)** jackstraw is a permute-the-scores shortcut and is likely anti-conservative, **(v)** the joint-rank `max(wedin, rand_dir, perm)` rule has no formal Type-I guarantee, and **(vi)** `assess_stability(target = "components")` ignores singular-vector sign/order ambiguity and will under-report stability. A concrete Verification Plan with runnable simulation snippets is included.
