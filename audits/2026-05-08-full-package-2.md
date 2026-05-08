# Scientific & Statistical Audit — `rajiveplus` (full package, second pass)

**Date:** 2026-05-08 (re-audit)
**Auditor mode:** read-only
**Scope:** primary analysis paths end-to-end. Independent re-audit; cross-references and confirms / updates the earlier report at [audits/2026-05-08-full-package.md](2026-05-08-full-package.md).
**Files inspected:** [R/Rajive.R](../R/Rajive.R), [R/Rajive_helpfunctions.R](../R/Rajive_helpfunctions.R), [R/RobustSVD.R](../R/RobustSVD.R), [R/jackstraw.R](../R/jackstraw.R), [R/simulation.functions.R](../R/simulation.functions.R), [R/visualization.R](../R/visualization.R), [src/RobustSVD.cpp](../src/RobustSVD.cpp), [src/RcppExports.cpp](../src/RcppExports.cpp), [tests/testthat/](../tests/testthat).

> **Read-only invariant honored.** No source / test / doc / config file was modified. No `devtools::document()`, install, or commit was run. Only this audit file was created (with `-2` suffix to avoid overwriting the existing same-day report).

---

## 1. Executive Summary

- **Overall verdict — partially defensible.** The core RaJIVE algebra (per-block robust SVD → stacked signal-score SVD → identifiability filter on per-block L2-projection magnitude → projector-based J / I / E split) is implemented coherently and matches the published algorithm (Ponzi et al. 2021; Feng et al. 2018). However, several quantitative claims (joint-rank thresholds, jackstraw p-values, bootstrap "stability" of components) rest on heuristic null distributions or shortcut estimators with no calibrated error-rate guarantee, and the parallel paths are not RNG-reproducible.
- **Resolved since prior audit:** prior Finding 2.1 (`sim_dist(num = 3)` undefined) is **fixed** — [R/simulation.functions.R:90-103](../R/simulation.functions.R#L90-L103) now implements `dist <- rexp(n*p, 1) - 1` plus an explicit `stop()` for unknown `num`. Five of the prior medium-risk items remain live; see §3.
- **Still high-risk (re-confirmed):**
  - `get_sv_threshold()` returns `NA` at the boundary `rank == length(d)` and `Rajive()`'s fast path triggers it deterministically when `initial_signal_ranks[k] >= min(dim(blocks[[k]]))` ([R/Rajive.R:115-118](../R/Rajive.R#L115-L118), [R/Rajive.R:67-86](../R/Rajive.R#L67-L86)).
  - Parallel RNG is not seeded with `L'Ecuyer-CMRG` / `doRNG`. With `num_cores > 1`, `set.seed()` does not reproduce Wedin / random-direction / permutation samples ([R/Rajive_helpfunctions.R:120-138, 152-189, 196-228](../R/Rajive_helpfunctions.R#L120-L228)).
- **New finding (high-risk):** `RobRSVD_all_cpp` re-uses the *original first-component* `(sinit1, uinit1, vinit1)` to **initialize every deflated rank-1 fit**, instead of the leading SVD of the residual `data − Red` ([src/RobustSVD.cpp:200-225](../src/RobustSVD.cpp#L200-L225)). For `nrank > 1` this gives a meaningless starting point for the M-estimator on the deflated matrix, which can stall on poor local solutions, slow convergence, or — when the residual's leading direction is nearly orthogonal to `uinit1` — converge to a non-leading direction. All downstream singular values / vectors at indices ≥ 2 inherit this bias.
- **New finding (medium-risk):** `RobRSVD1_cpp` freezes the scale estimate `mysigma = MAD(R₀) / 0.675` at iteration 0 and never updates it ([src/RobustSVD.cpp:54-59](../src/RobustSVD.cpp#L54-L59)). Standard Huber M-estimation re-estimates scale each iteration; the frozen-scale variant is biased when residuals shrink markedly during convergence.
- **New finding (medium-risk):** `showVarExplained_robust()` accesses singular values via positional `[[1]]` ([R/Rajive_helpfunctions.R:497-520](../R/Rajive_helpfunctions.R#L497-L520)), relying on the undocumented invariant that `RobRSVD_all_cpp` returns `list(d, u, v)` in that order. Swapping the C++ return order would silently break variance accounting (norms of `u`/`v` instead of `d`).
- **Still medium-risk (re-confirmed):** jackstraw is a permute-the-scores shortcut (anti-conservative); empirical p-value uses strict `>` not `≥`; `max(wedin, rand_dir, perm)` joint-rank rule has no calibrated FWER/FDR; Wedin samples bootstrap perp-basis columns rather than draw Haar-uniform; `assess_stability(target = "components")` ignores sign / order ambiguity; `associate_components(mode = "survival")` coerces `"alive"/"dead"` to `NA` via `as.logical()`.
- **Tests verify execution, not science.** No test in `tests/testthat/` checks (a) calibration of jackstraw p-values under the global null, (b) Type-I rate of joint-rank > 0 when no joint structure exists, or (c) recovery of true joint rank under known-signal simulation. The 207 reported passing tests confirm the API is internally consistent, not that the inferences are calibrated.

---

## 2. High-Risk Issues

### 2.1 `get_sv_threshold()` returns `NA` at the boundary `rank == length(d)`

- **Finding** — `get_sv_threshold(d, rank)` computes `0.5 * (d[rank] + d[rank+1])`. When `rank == length(d)`, `d[rank+1]` is `NA`, so the threshold is `NA`. The fast path in `Rajive()` constructs `svd_ranks <- pmin(initial_signal_ranks + 1L, max_rank_per_block)` and asks `get_svd_robustH(blocks[[k]], rank = svd_ranks[k])` for *exactly* that many singular values; whenever `initial_signal_ranks[k] + 1 > min(dim(blocks[[k]]))`, the returned `d` has length `initial_signal_ranks[k]`, then `get_sv_threshold(d, rank = initial_signal_ranks[k])` indexes `d[rank+1]` → `NA`.
- **Location** — definition at [R/Rajive.R:114-118](../R/Rajive.R#L114-L118); fast-path trigger at [R/Rajive.R:67-86](../R/Rajive.R#L67-L86); usage in identifiability filter at [R/Rajive.R:280-290](../R/Rajive.R#L280-L290) and individual-rank selection at [R/Rajive.R:325-340](../R/Rajive.R#L325-L340).
- **Evidence** — `singular_values <- mclapply(...)` then `sv_thresholds <- mapply(get_sv_threshold, singular_values, initial_signal_ranks)`. `NA` propagates: identifiability filter compares `if(sv < sv_thresholds[[k]])` → `if(NA)` → "missing value where TRUE/FALSE needed". Individual rank `sum(d > NA)` → `NA`, then `truncate_svd(rank = NA)` indexes `1:NA` → error.
- **Consequence** — `Rajive()` aborts (or yields NA-corrupted ranks) with an opaque error on small / square blocks, which is the standard regime when the user has supplied conservative `initial_signal_ranks` close to `min(dim)`. No user-facing diagnostic.

### 2.2 Parallel RNG is not configured

- **Finding** — `mclapply()` and `foreach::%dopar%` are used in the block SVD step ([R/Rajive.R:67, 76, 82](../R/Rajive.R#L67-L82)), Wedin resampling ([R/Rajive_helpfunctions.R:120-138](../R/Rajive_helpfunctions.R#L120-L138)), random-direction sampling ([R/Rajive_helpfunctions.R:165-189](../R/Rajive_helpfunctions.R#L165-L189)), and permutation bound sampling ([R/Rajive_helpfunctions.R:215-228](../R/Rajive_helpfunctions.R#L215-L228)). None of these set `RNGkind("L'Ecuyer-CMRG")`, `mc.set.seed = TRUE`, or use `doRNG::%dorng%`.
- **Evidence** — Worker processes inherit the parent RNG state at fork but advance independently. `rnorm(n_obs * l, ...)` inside random-direction `%dopar%` and `sample.int(n_obs)` inside permutation `%dopar%` will draw different sequences across runs even after `set.seed()` in the calling session.
- **Consequence** — Reproducibility — a load-bearing scientific property — silently fails for any user who runs `Rajive(..., num_cores > 1)`. Joint-rank thresholds and (in borderline cases) the joint rank itself can change run-to-run.

### 2.3 `RobRSVD_all_cpp` deflation re-uses original-component initialization for every rank

- **Finding** — Inside the `for (i in 1:(Rm-1))` deflation loop, every iteration calls `RobRSVD1_cpp(data - Red, sinit1, uinit1, vinit1, ...)` — using the **first** singular triplet of the **original** matrix as the warm start for the M-estimator on the deflated residual.
- **Location** — [src/RobustSVD.cpp:204-216](../src/RobustSVD.cpp#L204-L216):
  ```cpp
  for (int i = 1; i < Rm; i++) {
      List resi = RobRSVD1_cpp(data - Red, sinit1, uinit1, vinit1,
                                huberk, niter, tol);
      ...
  }
  ```
- **Consequence** —
  - The warm start `(uinit1, vinit1)` lies in the row/column space of the *first* component, which has just been subtracted out of `data - Red`. The M-estimator is therefore initialized at a point where the residual's leading direction is approximately *orthogonal* to the warm start.
  - In the best case this means many extra iterations to reach convergence; in the worst case the iteration converges to a non-leading direction or fails the `niter` cap silently (no convergence flag is returned).
  - All `d[2:nrank]`, `u[, 2:nrank]`, `v[, 2:nrank]` returned by `RobRSVD_all_cpp` are suspect for `nrank > 1`. This affects:
    - `get_svd_robustH(blocks[[k]], rank = signal_rank)` for `signal_rank > 1` — i.e. virtually every real use of `Rajive()`.
    - The stacked signal-score SVD `M_svd <- get_svd_robustH(M, rank = min(initial_signal_ranks))` when `min(initial_signal_ranks) > 1` — i.e. the matrix that *defines the joint scores*.
  - The R reference (`RobustSVD.R` in upstream RobRSVD packages) initializes deflated rank-1 fits from the leading `svd(data - Red)` triplet, not from the original `svdinit`. The C++ comment "matching the R source" is at odds with what a well-defined deflation strategy requires.
- **Disclosure note** — On uncontaminated data, `RobRSVD1_cpp` started from a poor warm start often still finds the right direction (the loss surface is benign); this likely explains why integration tests pass. Under contamination — the entire reason for using a robust SVD — the M-estimator is more sensitive to initialization, which is exactly when the fault matters most.

---

## 3. Medium-Risk Issues

### 3.1 `RobRSVD1_cpp` freezes the scale estimate `mysigma`

- **Finding** — `mysigma <- median(|vec(data - Appold)|) / 0.675` is computed once at iteration 0 from the initialization residual and never re-estimated inside the M-estimation loop.
- **Location** — [src/RobustSVD.cpp:54-59](../src/RobustSVD.cpp#L54-L59); the loop body re-uses the same `mysigma` for `Wmat = huberk / |R / mysigma|`.
- **Consequence** — As iterations progress, residuals typically shrink (well-fit case) or change distribution shape (under contamination). Standard Huber M-estimators re-estimate scale alongside location to keep the asymptotic efficiency claim. With frozen scale: when residuals shrink, every Huber weight saturates at 1 and the M-estimator silently degenerates to OLS on the deflated matrix — losing the robustness guarantee that motivates the package.
- **What would resolve the uncertainty** — A simulation under heavy-tailed (`dist.type = 3`) or mixture noise comparing recovered `u`, `v` against the R reference implementation that updates `mysigma` per iteration; or a controlled ε-contamination experiment where the breakdown of the frozen-scale variant should appear at lower ε than the iterated-scale version.

### 3.2 Jackstraw is a permute-the-scores shortcut

- **Finding** — Confirms prior audit §3.1. `generate_null_f_stats()` permutes the score vector against fixed feature rows rather than refitting Rajive after masking `s` synthetic null features. Disclosed in the docstring at [R/jackstraw.R:170-186](../R/jackstraw.R#L170-L186) but the calibration consequences are not.
- **Location** — [R/jackstraw.R:135-166](../R/jackstraw.R#L135-L166).
- **Evidence** — `null_f[, s] <- ols_f_stat_matrix(rows, scores_perm)` with `scores_perm <- sample(joint_comp_scores)`.
- **Consequence** — When a feature contributes appreciably to the estimated `joint_score`, permuting the score does not break the residual data-dependence between feature and score → null distribution is too narrow → p-values anti-conservative. Magnitude depends on per-feature loading × `1/n`.
- **What would resolve the uncertainty** — Calibration simulation; see Verification Plan 5.2.

### 3.3 Empirical p-value uses strict `>` rather than `≥`

- **Finding** — Confirms prior audit §3.2. `n_le <- findInterval(f_obs, sorted_pool)` (default `left.open = FALSE`) returns `#{j : pool[j] ≤ x}`, so `N - n_le = #{pool > x}`, not `#{pool ≥ x}` as the Phipson-Smyth formula and the docstring claim.
- **Location** — [R/jackstraw.R:97-110](../R/jackstraw.R#L97-L110).
- **Consequence** — Negligible for continuous F-statistics; non-negligible for small-`n` blocks with constant or near-constant features (more ties).
- **What would resolve the uncertainty** — Switch to `findInterval(f_obs, sorted_pool, left.open = TRUE)` (a one-line fix) and unit-test against `mean(pool >= f_obs)`.

### 3.4 Joint-rank `max(wedin_5%, rand_dir_95%, perm_95%)` rule has no calibrated error rate

- **Finding** — Confirms prior audit §3.3. The three thresholds have different statistical interpretations (one lower bound under H1, two upper bounds under H0). `max(..., na.rm = TRUE)` is then applied component-by-component to all squared singular values of `M_svd`, with no across-component multiple-testing correction.
- **Location** — [R/Rajive.R:201-231](../R/Rajive.R#L201-L231).
- **Consequence** — When `K` and per-block dimensions grow, the per-component Type-I rate for spurious joint components inflates; the joint rank is biased upward.
- **What would resolve the uncertainty** — Verification Plan 5.3 (Type-I simulation under H0).

### 3.5 Wedin bound bootstraps perp-basis columns instead of Haar-uniform sampling

- **Finding** — Confirms prior audit §3.4. `perp_resampled <- perp_basis[, sample.int(ncol, ncol, replace = TRUE)]` produces a non-orthonormal `r`-column matrix whose operator-2 norm is biased relative to the Haar-Wedin reference.
- **Location** — [R/Rajive_helpfunctions.R:120-138](../R/Rajive_helpfunctions.R#L120-L138).
- **Consequence** — Wedin samples are systematically larger than the proper Haar samples → the 5th-percentile threshold is too high → joint rank biased downward when Wedin dominates the `max()` rule.

### 3.6 `assess_stability(target = "components")` ignores sign / order ambiguity

- **Finding** — Confirms prior audit §3.5. The bootstrap correlation `cor(ref[idx, j], bs[, j])` is signed and assumes positional matching of bootstrap component `j` to reference component `j`; singular vectors carry arbitrary sign and are routinely permuted in repeat fits. Notably, the loadings branch *does* call `.procrustes_align()` ([R/visualization.R:1497-1499](../R/visualization.R#L1497-L1499)), so the asymmetry is internal inconsistency, not oversight in design.
- **Location** — [R/visualization.R:1551-1556](../R/visualization.R#L1551-L1556).
- **Consequence** — `mean_correlation` is biased toward 0; users see "stable" components reported as unstable.

### 3.7 `associate_components(mode = "survival")` silently coerces character status to NA

- **Finding** — Confirms prior audit §3.6. `survival::Surv(as.numeric(t_vec), as.logical(s_vec))` produces NA event indicators when `s_vec` is a character vector (`"alive"`, `"dead"`, etc.). No upstream validation.
- **Location** — [R/visualization.R:1320-1330](../R/visualization.R#L1320-L1330).
- **Consequence** — Cox model is fit on a mostly-censored Surv object; coefficients and p-values are meaningless without a stop / warning.

### 3.8 Variance accounting depends on undocumented positional invariant of robust-SVD return

- **Finding** — `showVarExplained_robust()` computes block-`k` joint variance as `norm(as.matrix(ajiveResults$block_decomps[[3*(k-1)+2]][[1]]), type = "F")^2 / norm(blocks[[k]], type = "F")^2`. The `[[1]]` accesses the first list element of the joint-decomposition list. This works *only* because `RobRSVD_all_cpp` happens to return `list(d = ..., u = ..., v = ...)` with `d` first, so `[[1]]` is the singular-value vector and `||d||_F^2 = sum(d^2) = ||J||_F^2`.
- **Location** — [R/Rajive_helpfunctions.R:497-520](../R/Rajive_helpfunctions.R#L497-L520); the underlying invariant is in [src/RobustSVD.cpp:218-223](../src/RobustSVD.cpp#L218-L223).
- **Consequence** — Reordering the C++ return list (e.g., to match base R's `svd()` convention `list(d, u, v)` becoming `list(u, d, v)` for any reason) would silently change `[[1]]` to `u` (an `n × rank` orthonormal matrix with `||u||_F^2 = rank`), making every reported variance fraction wrong without throwing an error.
- **What would resolve the uncertainty** — Replace `[[1]]` with `[["d"]]` (a one-line fix; not in scope for this read-only audit).

### 3.9 `truncate_svd(decomposition, rank = 0)` returns inconsistent shape

- **Finding** — When `rank == 0`, `truncate_svd()` writes back `u` as an `n × 1` zero matrix, `d` as scalar `0`, `v` as `p × 1` zero matrix ([R/Rajive_helpfunctions.R:39-50](../R/Rajive_helpfunctions.R#L39-L50)). Every other branch returns `n × rank` (zero-column when `rank == 0` would be the correct mirror).
- **Consequence** — Downstream code that does `ncol(decomp$u)` to determine the rank gets `1` instead of `0`. The companion zero-rank path in `get_joint_decomposition_robustH()` ([R/Rajive.R:325-360](../R/Rajive.R#L325-L360)) was patched to return zero-column matrices specifically to avoid this inconsistency, indicating the author already encountered the bug; the fix was not propagated to `truncate_svd()`.

### 3.10 `RobRSVD.all` evaluates `svd(data)` eagerly via a default argument

- **Finding** — `RobRSVD.all(data, nrank = min(dim(data)), svdinit = svd(data))` evaluates `svd(data)` whenever the caller does not pass `svdinit` — which is every internal call site. This is a full classical SVD on every robust-SVD invocation, including inside the random-direction `%dopar%` loop (now mitigated by switching that loop to base `svd()` directly — good) and inside every per-block call in `Rajive()`. Less of a correctness concern, more of a fragility issue: if `svd(data)` errors (e.g. non-finite entries), the user sees a confusing "error in svd(data) — infinite or missing values" before any `Rajive`-specific validation runs.
- **Location** — [R/RobustSVD.R:14](../R/RobustSVD.R#L14).
- **Consequence** — Weak; flagged for transparency. Edge-case error messages do not point at the data block responsible.

---

## 4. Unclear Assumptions Requiring Domain Review

1. **Block-wise centering / scaling expectation.** `Rajive()` does not center or rescale inputs, but the JIVE family universally assumes mean-zero blocks; the published RaJIVE paper (Ponzi et al. 2021, §2.1) further presumes Frobenius-norm rescaling so no single block dominates the stacked-signal SVD `M`. The vignettes do center but do not Frobenius-rescale. Should `Rajive()` (a) document this expectation prominently in `@details`, (b) perform centering by default with an `center = TRUE` argument, or (c) error if column means are not near zero?
2. **Frozen-scale Huber (Finding 3.1).** Is the deviation from iterated-scale Huber an intentional efficiency trade-off, or a subtle bug? The user-facing docs do not mention the choice.
3. **Permutation-bound exchangeability (carried over).** `get_perm_bound_robustH()` independently permutes block rows. For paired / longitudinal designs this is the wrong null. Should the API accept a `strata` argument?
4. **Jackstraw target of inference.** The docstring at [R/jackstraw.R:185-189](../R/jackstraw.R#L185-L189) correctly notes the test is on the *estimated* score, not the latent one. Is that distinction acceptable for the manuscript's stated use cases (microbiome, CLL multi-omics)?
5. **BH multiple-testing across highly-correlated features.** All multiple-testing in `jackstraw_rajive()` defaults to `BH`, which is conservative under positive regression dependency but can be invalid under arbitrary dependency. For omics features with strong block-internal correlation, is `BY` the more defensible default? At minimum, document the assumption.

---

## 5. Verification Plan

For each finding, a runnable check (do not modify the package; run from a scratch script in `R4_51`).

### 5.1 Boundary failure of `get_sv_threshold()` (Finding 2.1)

```r
library(rajiveplus)
set.seed(1)
X   <- matrix(rnorm(20), 5, 4)              # min(dim) = 4
fit <- rajiveplus:::get_svd_robustH(X, rank = 4)
th  <- rajiveplus:::get_sv_threshold(fit$d, rank = 4)
testthat::expect_false(is.na(th))           # expected FAIL today (returns NA)

# End-to-end trigger:
Y <- ajive.data.sim(K = 2, rankJ = 1, rankA = c(2, 2), n = 5, pks = c(4, 4))
testthat::expect_no_error(
  Rajive(Y$sim_data, initial_signal_ranks = c(4, 4), joint_rank = 1)
)                                           # expected FAIL today
```

### 5.2 Parallel RNG reproducibility (Finding 2.2)

```r
Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = 40, pks = c(30, 20))
set.seed(7); a <- Rajive(Y$sim_data, c(5,4), n_rand_dir_samples = 100, num_cores = 2)
set.seed(7); b <- Rajive(Y$sim_data, c(5,4), n_rand_dir_samples = 100, num_cores = 2)
testthat::expect_equal(a$joint_rank_sel$rand_dir$rand_dir_samples,
                       b$joint_rank_sel$rand_dir$rand_dir_samples)
# Expected FAIL today; PASS after L'Ecuyer-CMRG / doRNG wiring.
```

### 5.3 RobRSVD deflation initialization (Finding 2.3)

```r
# Property: with nrank = 3 on a clean rank-3 matrix, RobRSVD_all_cpp should
# recover singular values matching base svd() to high precision.
set.seed(42)
U  <- qr.Q(qr(matrix(rnorm(60), 20, 3)))
V  <- qr.Q(qr(matrix(rnorm(60), 20, 3)))
D  <- diag(c(10, 5, 1))
X  <- U %*% D %*% t(V) + matrix(rnorm(400, sd = 0.01), 20, 20)

ref <- svd(X, nu = 3, nv = 3)$d[1:3]
rob <- rajiveplus:::RobRSVD.all(X, nrank = 3)$d
abs(rob - ref) / ref
# If components 2 and 3 deviate by >> 1% (with this SNR), the deflation
# warm-start bug is biting.  Compare also the angle between rob$u[,2:3]
# and ref$u[,2:3] via principal angles (subspace::angle()).
```

### 5.4 Frozen-scale Huber sensitivity (Finding 3.1)

```r
# Inject row-wise outliers and compare RobRSVD1_cpp leading direction
# against an iterated-scale reference (e.g., MASS::rlm or pcaPP::PCAproj).
set.seed(1)
X <- matrix(rnorm(200), 20, 10); X[1, ] <- X[1, ] + 50
init <- svd(X, nu = 1, nv = 1)
rob  <- rajiveplus:::RobRSVD1(X, sinit = init$d[1],
                              uinit = init$u[, 1], vinit = init$v[, 1])
ref  <- pcaPP::PCAproj(X, k = 1, scale = NULL)$loadings[, 1]
abs(crossprod(rob$v, ref))                  # should be near 1
```

### 5.5 Calibration of jackstraw p-values under H0 (Finding 3.2)

```r
set.seed(2026); B <- 200
ps <- replicate(B, {
  blk <- list(matrix(rnorm(60*40), 60, 40), matrix(rnorm(60*30), 60, 30))
  fit <- Rajive(blk, c(5, 4), joint_rank = 1, num_cores = 1L)
  js  <- jackstraw_rajive(fit, blk, n_null = 50, correction = "none")
  unlist(lapply(js, function(b) lapply(b, `[[`, "p_values")))
})
ks.test(ps, "punif")$p.value      # expect > 0.1 under calibration
mean(ps <= 0.05)                   # expect ~ 0.05; >> 0.05 confirms anti-conservativeness
```

### 5.6 Empirical p-value `>=` vs `>` (Finding 3.3)

```r
f_obs  <- c(1, 2, 3)
f_null <- matrix(c(1, 1, 2, 2, 3, 3), nrow = 3)
expected <- (1 + c(6, 4, 2)) / (1 + 6)        # using >=
got      <- rajiveplus:::compute_empirical_pvalues(f_obs, f_null)
testthat::expect_equal(got, expected)         # expected FAIL today
```

### 5.7 Joint-rank Type-I rate under H0 (Finding 3.4)

```r
set.seed(2026)
ranks <- replicate(200, {
  blk <- lapply(c(50, 40, 30), function(p) matrix(rnorm(60*p), 60, p))
  Rajive(blk, c(6, 5, 4), n_wedin_samples = 200,
         n_rand_dir_samples = 200, num_cores = 1L)$joint_rank
})
mean(ranks > 0)            # nominal 0.05; report observed
table(ranks)
# Repeat with n_perm_samples = 200 in place of n_rand_dir_samples.
```

### 5.8 Wedin Haar-vs-bootstrap discrepancy (Finding 3.5)

```r
# KS-test the current bootstrap-of-columns wedin samples against a
# Haar-uniform reference: Q ~ qr.Q(qr(matrix(rnorm(p*r), p, r))) and
# project through perp_basis.  Large D confirms bias.
```

### 5.9 `assess_stability(components)` sign / order invariance (Finding 3.6)

```r
set.seed(7)
Y    <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = 60, pks = c(40, 30))
fit  <- Rajive(Y$sim_data, c(5, 4))
stab <- assess_stability(fit, Y$sim_data, c(5, 4),
                         target = "components", B = 20L, sample_frac = 0.9)
# Today: many entries near 0 due to sign flips.
# Property check after fix: all(abs(stab$mean_correlation) > 0.7)
```

### 5.10 Survival status coercion (Finding 3.7)

```r
meta <- data.frame(t = runif(40, 1, 10),
                   status = sample(c("alive","dead"), 40, TRUE))
fit  <- Rajive(ajive.data.sim(K=2, rankJ=2, rankA=c(5,4), n=40, pks=c(30,20))$sim_data,
               c(5,4))
res  <- associate_components(fit, meta, mode = "survival",
                             time_col = "t", status_col = "status",
                             variable = "t")
# Today: silently fits an all-censored Cox model; res$p_value is meaningless.
```

### 5.11 Variance-accounting positional invariant (Finding 3.8)

```r
fit <- Rajive(ajive.data.sim(K=2, rankJ=2, rankA=c(5,4), n=40, pks=c(30,20))$sim_data,
              c(5,4))
ve  <- showVarExplained_robust(fit, blocks)
testthat::expect_true(all(unlist(ve$Joint) >= 0 & unlist(ve$Joint) <= 1))
testthat::expect_true(all(unlist(ve$Resid) >= -1e-8))   # exposes negative residual var
# Then spot-check that VarJoint matches recomputation via [["d"]] not [[1]]:
J_d <- fit$block_decomps[[2]][["d"]]
sum(J_d^2) / norm(blocks[[1]], type = "F")^2     # should equal ve$Joint[[1]]
```

### 5.12 Recovery under known signal (sanity)

```r
set.seed(2026)
hits <- replicate(50, {
  Y <- ajive.data.sim(K = 3, rankJ = 2, rankA = c(5, 4, 3),
                      n = 80, pks = c(60, 50, 40), dist.type = 1)
  Rajive(Y$sim_data, c(5, 4, 3))$joint_rank
})
mean(hits == 2)        # expect close to 1 if the pipeline works on clean data
```

---

## 6. Verification Standard / Trace

### 6.1 Trace: `Rajive()` joint-rank selection path

```
user → Rajive(blocks, initial_signal_ranks,
              n_wedin_samples=1000, n_rand_dir_samples=1000,
              joint_rank=NA, n_perm_samples=NA, num_cores=1L)        R/Rajive.R:50
  ├─ block SVD: get_svd_robustH(X_k)  (mclapply)                     R/Rajive.R:67/82
  │     └─ RobRSVD.all → RobRSVD_all_cpp                              R/RobustSVD.R:14
  │         └─ deflation reuses (sinit1,uinit1,vinit1)  ⚠ Finding 2.3 src/RobustSVD.cpp:204
  │             └─ RobRSVD1_cpp: frozen mysigma          ⚠ Finding 3.1 src/RobustSVD.cpp:54
  ├─ get_sv_threshold(d_k, rank=initial_signal_ranks[k])              R/Rajive.R:90
  │       ⚠ NA at boundary (Finding 2.1)                              R/Rajive.R:114
  ├─ get_joint_scores_robustH(...)                                    R/Rajive.R:97
  │   ├─ stack signal scores → M; M_svd                                R/Rajive_helpfunctions(.R+inline)
  │   ├─ get_wedin_bound_samples (per block)                           R/Rajive_helpfunctions.R:88
  │   │     └─ wedin_bound_resampling  ⚠ bootstrap not Haar (Finding 3.5)
  │   │         ⚠ parallel RNG (Finding 2.2)
  │   ├─ get_perm_bound_robustH  OR  get_random_direction_bound_robustH
  │   │     ⚠ parallel RNG (Finding 2.2)
  │   ├─ overall_t = max(wedin_5%, rand_95%, perm_95%, na.rm=TRUE)
  │   │     ⚠ no calibrated FWER/FDR (Finding 3.4)
  │   └─ identifiability filter (L2-norm of t(blocks[[k]]) %*% u_j)    R/Rajive.R:265-285
  └─ get_final_decomposition_robustH (per block, mapply)               R/Rajive.R:111
        ├─ get_individual_decomposition_robustH   (uses sv_thresholds, possibly NA)
        └─ get_joint_decomposition_robustH        (zero-rank guard at L:330-345)
```

### 6.2 Trace: `jackstraw_rajive()`

```
user → jackstraw_rajive(ajive_output, blocks, alpha, n_null, correction)  R/jackstraw.R:235
  └─ for k in 1..K:
       for j in 1..joint_rank:
         f_obs  ← ols_f_stat_matrix(X_t, scores_j)                         R/jackstraw.R:299
         f_null ← generate_null_f_stats(X_t, scores_j, n_null)             R/jackstraw.R:300
                  ⚠ permute scores, not refit (Finding 3.2)
         p_vals ← compute_empirical_pvalues(f_obs, f_null)                 R/jackstraw.R:301
                  ⚠ '>' not '>=' (Finding 3.3)
  └─ all_p_adj ← p.adjust(unlisted, method = correction)                   R/jackstraw.R:329-336
```

### 6.3 Claim provenance

- **Test totals (207 passing per [PROGRESS.md](../PROGRESS.md))** verify *execution* and *consistency with previous behaviour*. They do **not** verify statistical calibration. None of `tests/testthat/test-jackstraw.R`, `test-rajive-perm.R`, or `test-rajive-rcpp.R` runs a null simulation; the calibration plan in §5.5–§5.7 is essential and not currently exercised.
- **Robustness claims** (the package's reason for existing) are not verified anywhere: there is no test that contaminates inputs and confirms the recovered subspace is closer to truth than what classical aJIVE would give. Until §5.4 is run, the robustness claim rests on documentation alone.

---

## Pointer / Executive Summary for chat

Audit report (second pass) written to [audits/2026-05-08-full-package-2.md](audits/2026-05-08-full-package-2.md).

Headline: the prior audit's `sim_dist(num=3)` finding is **fixed**; everything else still holds. **New high-risk:** `RobRSVD_all_cpp` deflation re-uses the *original* matrix's leading singular triplet to warm-start the M-estimator on every deflated rank, which is the wrong starting point and can corrupt all singular values / vectors at indices ≥ 2 — the regime that defines the joint scores for any `min(initial_signal_ranks) > 1`. **New medium-risk:** `RobRSVD1_cpp` freezes the scale estimate `mysigma` at iteration 0 (degrading Huber robustness when residuals shrink); `showVarExplained_robust()` relies on undocumented positional `[[1]]` access into the robust-SVD return list; `truncate_svd(rank = 0)` returns a 1-column-zero matrix instead of a 0-column matrix. Carried over: `get_sv_threshold()` boundary NA, parallel RNG not seeded, jackstraw permute-the-scores shortcut (anti-conservative), empirical p-value `>` vs `≥`, `max(wedin, rand_dir, perm)` rule with no calibrated error rate, Wedin bootstrap-of-columns not Haar, `assess_stability(components)` ignores sign / order, survival status coercion silent-NA. A runnable Verification Plan is included for every finding.
