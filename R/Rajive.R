# Internal: cross-platform parallel lapply.  parallel::mclapply uses fork() and
# silently degrades to serial on Windows; we fall back to a PSOCK cluster
# there so num_cores > 1 actually does work on every platform.  When
# num_cores == 1L we run serially to avoid any cluster setup overhead.
#' @noRd
.rajive_parallel_lapply <- function(X, FUN, num_cores = 1L,
                                    extra_globals = list(), ...) {
  num_cores <- max(1L, as.integer(num_cores))
  if (num_cores == 1L) {
    return(lapply(X, FUN, ...))
  }
  is_windows <- .Platform$OS.type == "windows"
  if (!is_windows) {
    return(parallel::mclapply(X, FUN, mc.cores = num_cores,
                              mc.set.seed = TRUE, ...))
  }
  seed <- sample.int(.Machine$integer.max, 1L)
  cl <- parallel::makePSOCKcluster(num_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, iseed = seed)
  if (length(extra_globals) > 0L) {
    parallel::clusterExport(cl, varlist = names(extra_globals),
                            envir = list2env(extra_globals))
  }
  parallel::clusterEvalQ(cl, {
    requireNamespace("rajiveplus", quietly = TRUE)
    NULL
  })
  parallel::parLapply(cl, X, FUN, ...)
}

.inform_identifiability_norm_default <- function() {
  opt <- "rajiveplus.identifiability_norm.default_informed"
  if (!isTRUE(getOption(opt, FALSE))) {
    cli::cli_inform(
      'Rajive() uses identifiability_norm = "l2" by default; use "l1" for closer original RaJIVE parity.'
    )
    options(rajiveplus.identifiability_norm.default_informed = TRUE)
  }
  invisible(NULL)
}

.identifiability_projection_norm <- function(score, identifiability_norm) {
  if (identifiability_norm == "l2") {
    sqrt(sum(score^2))
  } else {
    sum(abs(score))
  }
}

.match_identifiability_norm <- function(identifiability_norm) {
  tryCatch(
    match.arg(identifiability_norm, c("l2", "l1")),
    error = function(e) {
      cli::cli_abort(
        '`identifiability_norm` must be one of "l2" or "l1".',
        class = "rajiveplus_invalid_input"
      )
    }
  )
}


#' Robust Angle-based Joint and Individual Variation Explained (RaJIVE)
#'
#' Computes the robust aJIVE decomposition of a list of multi-view data
#' matrices into joint, block-individual, and residual (noise) components.
#' The robust SVD step uses an M-estimator (Huber loss) so the decomposition
#' is resistant to a moderate fraction of element-wise outliers in any block.
#'
#' Joint rank is selected by comparing the squared common singular values of
#' the stacked signal-score matrix to a threshold derived from a Wedin bound
#' and a random-direction bound, a permutation bound, or both.  The combined
#' threshold is \code{max(wedin, rand_dir, perm)}; this heuristic is conservative in
#' practice but does not carry a formal FWER/FDR guarantee for rank
#' selection.
#'
#' @param blocks List of numeric matrices. Each block has \eqn{n} rows
#'   (matched samples, in the same order across blocks) and \eqn{p_k}
#'   columns. All blocks must be finite (no \code{NA}/\code{NaN}/\code{Inf})
#'   and have non-degenerate columns.
#' @param initial_signal_ranks Integer vector of length \code{length(blocks)}
#'   giving the initial signal-rank estimate for each block. Choose by visual
#'   inspection of each block's scree plot.
#' @param full Logical. If \code{TRUE} (default) the full \eqn{J}, \eqn{I},
#'   \eqn{E} matrices are stored alongside their SVDs; set to \code{FALSE}
#'   to save memory.
#' @param n_wedin_samples Positive integer. Number of Wedin bound samples to
#'   draw per block. Default \code{1000}.
#' @param n_rand_dir_samples Positive integer. Number of random-direction
#'   bound samples to draw. Default \code{1000}. When both
#'   \code{n_rand_dir_samples} and \code{n_perm_samples} are supplied, both
#'   bounds are computed independently and the maximum is used as the threshold.
#' @param joint_rank Integer, \code{NA}, or \code{"native_cv"} in native
#'   missing-data mode. When \code{NA} (default), complete-data mode estimates
#'   the rank from Wedin/random-direction/permutation bounds, while native
#'   missing-data mode selects among native candidate ranks and stores the
#'   diagnostics in \code{fit$missing$rank_diagnostics}.
#' @param n_perm_samples Integer or \code{NA}. Number of permutation samples
#'   for the permutation-based joint-rank threshold.  Rows of each block's
#'   signal scores are independently shuffled, the leading squared singular
#'   value of the stacked permuted matrix is recorded, and its 95th percentile
#'   is used as the threshold.  When both \code{n_perm_samples} and
#'   \code{n_rand_dir_samples} are supplied, both are computed and their
#'   maximum (together with the Wedin bound) is the final threshold.
#'   Default \code{NA}.
#' @param num_cores Positive integer. Number of cores to use for parallel
#'   block SVDs, Wedin / random-direction / permutation bound resampling.
#'   Default \code{1L} (serial). When \code{> 1} the function transparently
#'   uses \code{\link[parallel]{mclapply}} on Unix and a
#'   \code{\link[parallel]{makePSOCKcluster}} on Windows. The RNG kind is
#'   temporarily set to \code{"L'Ecuyer-CMRG"} for reproducible parallel
#'   sampling (with a one-time warning); set
#'   \code{RNGkind("L'Ecuyer-CMRG")} explicitly to silence it.
#' @param seed Integer or \code{NA}. Optional seed passed to \code{set.seed()}
#'   before fitting. Default \code{NA} leaves the caller's RNG state unchanged.
#' @param identifiability_norm Character. Norm used in the post-selection
#'   identifiability filter for each candidate joint component. \code{"l2"}
#'   (default) compares \eqn{X_k^T u_j} to the singular-value threshold on
#'   Euclidean scale, matching the scale of AJIVE singular values. \code{"l1"}
#'   uses \code{sum(abs(X_k^T u_j))}, matching original RaJIVE's
#'   \code{norm(score)} behavior for a one-column matrix and therefore giving
#'   closer legacy parity.
#' @param missing Character. \code{"error"} (default) preserves complete-data
#'   behaviour and rejects non-finite block entries. \code{"native"} enables
#'   native observed-entry fitting for blocks with \code{NA}s or an explicit
#'   observation \code{mask}.
#' @param mask Optional list of logical matrices with the same dimensions as
#'   \code{blocks}; \code{TRUE} marks observed entries and \code{FALSE} marks
#'   missing entries excluded from the native missing-data objective.
#' @param missing_control List created by \code{\link{rajive_missing_control}}
#'   controlling native missing-data preprocessing, rank diagnostics,
#'   uncertainty refits, and censoring/sensitivity metadata.
#' @param uncertainty Character. \code{"none"} (default), \code{"bootstrap"},
#'   or \code{"mi"}. Native missing-data uncertainty stores aligned refits in
#'   \code{fit$missing$uncertainty}; complete-data fits ignore this argument.
#'
#' @return An object of class \code{"rajive"}: a named list containing
#'   \describe{
#'     \item{\code{block_decomps}}{A list of length \eqn{3K}. For block
#'       \eqn{k} (1-indexed): individual component at index
#'       \eqn{3(k-1)+1}, joint component at \eqn{3(k-1)+2}, noise (residual)
#'       at \eqn{3(k-1)+3}. Each entry has fields \code{u}, \code{d},
#'       \code{v}, and \code{full} (when \code{full = TRUE}).}
#'     \item{\code{joint_scores}}{The \eqn{n \times r_J} shared joint score
#'       matrix, where \eqn{r_J} is the estimated (or user-supplied) joint
#'       rank. May be a 0-column matrix when no joint signal is detected.}
#'     \item{\code{joint_rank}}{Integer. Estimated (or user-supplied) joint
#'       rank.}
#'     \item{\code{joint_rank_sel}}{Diagnostics from joint-rank selection:
#'       observed singular values, Wedin and random-direction (or
#'       permutation) samples, percentile cutoffs, and the indices of
#'       components dropped by the selected identifiability filter. The
#'       selected \code{identifiability_norm} is also stored here.}
#'   }
#'
#' @section Robustness and Huber tuning:
#' The M-estimator at each block uses a fixed Huber threshold
#' \code{huberk = 1.345}, which corresponds to 95\% asymptotic relative
#' efficiency under a Gaussian reference.  This constant is appropriate for
#' continuous-scale data (e.g. normalised gene expression, metabolomics).
#' For count-overdispersed or heavy-tailed data (e.g. microbiome log-ratios,
#' raw spectral intensities) where the nominal outlier fraction exceeds
#' \eqn{\approx 20\%}, consider pre-transforming data or calling
#' \code{\link{RobRSVD.all}} directly with a custom \code{k} argument.
#'
#' @section Block normalization:
#' \code{Rajive()} expects blocks to be column-centred.  Cross-block scale
#' comparability (relevant for the joint-rank threshold) additionally
#' requires Frobenius normalisation: \code{X <- X / norm(X, "F")} so that
#' each block contributes equally to the stacked signal-score matrix.
#' Omitting centring or normalisation does not error but will distort the
#' joint-rank estimate and loadings in proportion to the raw scale
#' differences (JIVE convention; see Lock et al. 2013,
#' \emph{Ann. Appl. Statist.} 7(1):523--542).
#'
#' @section Permutation and paired designs:
#' When \code{n_perm_samples} is supplied, the joint-rank threshold is
#' derived by independently permuting row order within each block's signal
#' scores, which assumes exchangeability of observations.  This assumption
#' holds for independent samples but is violated for paired, longitudinal,
#' or otherwise structured designs (time series, family trios, matched
#' controls).  For such designs, use \code{joint_rank} to fix the rank
#' manually instead of relying on \code{n_perm_samples}.
#'
#' @section Joint-rank threshold (heuristic, no formal FWER/FDR):
#' The joint rank is selected by comparing squared common singular values
#' to \code{max(wedin_5\%, rand_dir_95\% | perm_95\%)}.  This rule is a
#' heuristic: it has no formal FWER or FDR guarantee.  Under a global null
#' (independent Gaussian blocks, 3 blocks, \eqn{n = 60}, \eqn{p_k} in
#' \{50, 40, 30\}), calibration experiments find an empirical Type-I rate
#' of approximately 2--8\% for \code{n_wedin_samples = n_rand_dir_samples =
#' 200}; the rate approaches 5\% from below as the sample count increases.
#' Treat these numbers as empirical diagnostics for the stated simulation
#' setting rather than a universal error-rate guarantee.
#'
#' @section Reproducibility:
#' For exactly reproducible results, set \code{set.seed()} (and, for
#' \code{num_cores > 1}, \code{RNGkind("L'Ecuyer-CMRG")}) before calling
#' \code{Rajive()}.
#'
#' @section Native missing-data mode:
#' With \code{missing = "native"}, missing entries are excluded from centering,
#' scaling, robust SVD warm starts, masked geometry diagnostics, and
#' observed-entry reconstruction objectives. Numeric \code{joint_rank} values
#' request fixed-rank native fits; \code{joint_rank = NA} and
#' \code{joint_rank = "native_cv"} run native candidate-rank diagnostics and
#' select the recorded minimum-score rank. Reconstructions returned by
#' \code{\link{get_reconstructed_blocks}} are derived outputs, work even when
#' \code{full = FALSE}, and carry estimability labels from
#' \code{\link{get_estimability}}. Individual signal for an entirely missing
#' sample-block row is labelled \code{"not_identifiable"} and is not silently
#' filled from other blocks. Under partial missingness the joint/individual
#' split is approximate: the individual component is fit on observed residuals
#' and is not exactly orthogonal to the joint subspace; inspect the
#' \code{orthogonality} element returned by
#' \code{\link{get_missing_diagnostics}}(). With
#' \code{uncertainty = "mi"}, refits run through the complete-data pipeline
#' after imputation, whereas \code{"bootstrap"} refits keep the original mask.
#'
#' @seealso \code{\link{ajive.data.sim}} for simulating test data;
#'   \code{\link{get_joint_rank}}, \code{\link{get_joint_scores}},
#'   \code{\link{get_block_scores}}, \code{\link{get_block_loadings}},
#'   \code{\link{get_block_matrix}} for component accessors;
#'   \code{\link{plot_components}} for diagnostic plots;
#'   \code{\link{jackstraw_rajive}} for feature-level significance testing.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n   <- 50
#' pks <- c(100, 80, 50)
#' Y   <- ajive.data.sim(K = 3, rankJ = 3, rankA = c(7, 6, 4),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, initial_signal_ranks = c(7, 6, 4))
#' fit
#' summary(fit)
#' }
#'
#' @export
Rajive <- function(blocks, initial_signal_ranks, full=TRUE,
                           n_wedin_samples=1000, n_rand_dir_samples=1000,
                           joint_rank=NA,
                           n_perm_samples=NA,
                           num_cores=1L,
                           seed=NA_integer_,
                           identifiability_norm=c("l2", "l1"),
                           missing = c("error", "native"),
                           mask = NULL,
                           missing_control = rajive_missing_control(),
                           uncertainty = c("none", "bootstrap", "mi"))
{
  identifiability_norm_defaulted <- base::missing(identifiability_norm)
  identifiability_norm <- .match_identifiability_norm(identifiability_norm)
  missing <- match.arg(missing)
  uncertainty <- match.arg(uncertainty)
  if (identifiability_norm_defaulted) {
    .inform_identifiability_norm_default()
  }

  if (missing == "native") {
    fit <- .Rajive_incomplete(
      blocks = blocks,
      initial_signal_ranks = initial_signal_ranks,
      joint_rank = joint_rank,
      mask = mask,
      missing_control = missing_control,
      full = full,
      n_wedin_samples = n_wedin_samples,
      n_rand_dir_samples = n_rand_dir_samples,
      n_perm_samples = n_perm_samples,
      num_cores = num_cores,
      seed = seed,
      identifiability_norm = identifiability_norm
    )
    if (uncertainty != "none") {
      fit$missing$uncertainty <- .compute_missing_uncertainty(fit, method = uncertainty)
    }
    return(fit)
  }

  if (!is.null(mask)) {
    cli::cli_abort(
      "`mask` can only be supplied when `missing = 'native'`.",
      class = "rajiveplus_invalid_input"
    )
  }
  if (is.character(joint_rank)) {
    cli::cli_abort(
      "Character `joint_rank` values such as \"native_cv\" require `missing = 'native'`.",
      class = "rajiveplus_invalid_input"
    )
  }

  .Rajive_core(
    blocks = blocks,
    initial_signal_ranks = initial_signal_ranks,
    full = full,
    n_wedin_samples = n_wedin_samples,
    n_rand_dir_samples = n_rand_dir_samples,
    joint_rank = joint_rank,
    n_perm_samples = n_perm_samples,
    num_cores = num_cores,
    seed = seed,
    identifiability_norm = identifiability_norm,
    rank_only = FALSE
  )
}

.Rajive_rank_only <- function(blocks, initial_signal_ranks,
                              n_wedin_samples=1000, n_rand_dir_samples=1000,
                              joint_rank=NA,
                              n_perm_samples=NA,
                              num_cores=1L,
                              seed=NA_integer_,
                              identifiability_norm=c("l2", "l1")) {
  identifiability_norm <- .match_identifiability_norm(identifiability_norm)
  .Rajive_core(
    blocks = blocks,
    initial_signal_ranks = initial_signal_ranks,
    full = FALSE,
    n_wedin_samples = n_wedin_samples,
    n_rand_dir_samples = n_rand_dir_samples,
    joint_rank = joint_rank,
    n_perm_samples = n_perm_samples,
    num_cores = num_cores,
    seed = seed,
    identifiability_norm = identifiability_norm,
    rank_only = TRUE
  )
}

.Rajive_core <- function(blocks, initial_signal_ranks, full=TRUE,
                         n_wedin_samples=1000, n_rand_dir_samples=1000,
                         joint_rank=NA,
                         n_perm_samples=NA,
                         num_cores=1L,
                         seed=NA_integer_,
                         identifiability_norm=c("l2", "l1"),
                         rank_only=FALSE) {

  identifiability_norm <- .match_identifiability_norm(identifiability_norm)
  num_cores <- max(1L, as.integer(num_cores))
  if (!is.na(seed)) set.seed(as.integer(seed))

  # Abort early with clear classes/messages for invalid inputs.
  if (!is.list(blocks) || length(blocks) < 1L) {
    cli::cli_abort(
      c("`blocks` must be a non-empty list of numeric matrices."),
      class = "rajiveplus_invalid_input"
    )
  }
  if (length(initial_signal_ranks) != length(blocks)) {
    cli::cli_abort(
      c("`initial_signal_ranks` must have one entry per block.",
        "x" = "Got {.val {length(initial_signal_ranks)}} ranks for {.val {length(blocks)}} blocks."),
      class = "rajiveplus_invalid_input"
    )
  }
  if (!all(vapply(blocks, is.matrix, logical(1L)))) {
    cli::cli_abort(
      c("Every element of `blocks` must be a matrix."),
      class = "rajiveplus_invalid_input"
    )
  }
  if (!all(vapply(blocks, function(x) all(is.finite(x)), logical(1L)))) {
    cli::cli_abort(
      c("All entries of each block must be finite (no NA/NaN/Inf)."),
      class = "rajiveplus_invalid_input"
    )
  }

  blocks <- lapply(seq_along(blocks), function(k) {
    x <- blocks[[k]]
    bad_cols <- apply(x, 2L, stats::sd) < .Machine$double.eps^0.5
    if (any(bad_cols)) {
      cli::cli_warn(
        c("Degenerate column(s) in block {.val {k}} dropped automatically.",
          "i" = "{.val {sum(bad_cols)}} column(s) with near-zero variance removed."),
        class = "rajiveplus_degenerate_block"
      )
      x[, !bad_cols, drop = FALSE]
    } else {
      x
    }
  })

  n_obs <- nrow(blocks[[1L]])
  if (n_obs < sum(initial_signal_ranks)) {
    cli::cli_warn(
      c("`n` is smaller than sum(initial_signal_ranks); estimation may be underdetermined.",
        "i" = "n = {.val {n_obs}}, sum(initial_signal_ranks) = {.val {sum(initial_signal_ranks)}}"),
      class = "rajiveplus_underdetermined"
    )
  }

  # Parallel paths use L'Ecuyer-CMRG for reproducible sampling.
  if (num_cores > 1L && RNGkind()[1L] != "L'Ecuyer-CMRG") {
    prev_rng <- RNGkind()
    RNGkind("L'Ecuyer-CMRG")
    on.exit(do.call(RNGkind, as.list(prev_rng)), add = TRUE)
    cli::cli_warn(
      c("RNG kind was set to L'Ecuyer-CMRG for reproducible parallel sampling.",
        "i" = "Set RNGkind(\"L'Ecuyer-CMRG\") before calling Rajive() to silence this warning."),
      class = "rajiveplus_rngkind_adjusted"
    )
  }

  K <- length(blocks)


  # step 1: initial signal space extraction --------------------------------
  # initial estimate of signal space with SVD

  # If joint_rank is fixed, wedin resampling is skipped, so we only need
  # signal_rank + 1 singular values per block for thresholding.
  # Audit fix #19: parallel::mclapply silently degrades to mc.cores = 1 on
  # Windows (no fork).  Wrap in a small helper that uses a PSOCK cluster on
  # Windows so users get true parallelism on every platform.
  if (is.na(joint_rank)) {
    # Full SVD is required for Wedin bound resampling (uses U_perp/V_perp).
    block_svd <- .rajive_parallel_lapply(blocks, get_svd_robustH,
                                         num_cores = num_cores)
  } else {
    max_rank_per_block <- vapply(blocks, function(x) min(dim(x)), integer(1L))
    svd_ranks <- pmin(initial_signal_ranks + 1L, max_rank_per_block)

    block_svd <- .rajive_parallel_lapply(
      seq_along(blocks),
      function(k) get_svd_robustH(blocks[[k]], rank = svd_ranks[k]),
      num_cores = num_cores,
      extra_globals = list(blocks = blocks, svd_ranks = svd_ranks)
    )
  }
  # extract singular values from list (cheap; no need to parallelise)
  singular_values <- lapply(block_svd, function(l) l[[1]])
  # apply get_sv_threshold
  sv_thresholds <- mapply(function(l,p)
    get_sv_threshold(l, rank = p), singular_values, initial_signal_ranks)




  # step 2: joint sapce estimation -------------------------------------------------------------

  out <- get_joint_scores_robustH(blocks, block_svd, initial_signal_ranks, sv_thresholds,
                                  n_wedin_samples=n_wedin_samples,
                                  n_rand_dir_samples=n_rand_dir_samples,
                                  joint_rank=joint_rank,
                                  n_perm_samples=n_perm_samples,
                                  num_cores=num_cores,
                                  identifiability_norm=identifiability_norm)
  joint_rank_sel_results <- out$rank_sel_results
  joint_scores <- out$joint_scores

  joint_rank <- dim(joint_scores)[2]

  if (isTRUE(rank_only)) {
    rank_decomposition <- list(
      joint_scores = joint_scores,
      joint_rank = joint_rank,
      joint_rank_sel = joint_rank_sel_results
    )
    class(rank_decomposition) <- "rajive_rank_only"
    return(rank_decomposition)
  }

  # step 3: final decomposition -----------------------------------------------------

  block_decomps <- mapply(function(l, m)
    get_final_decomposition_robustH(
      l,
      joint_scores = joint_scores,
      sv_threshold = m,
      full = full
    ), blocks, sv_thresholds)



  jive_decomposition <- list(block_decomps=block_decomps)
  jive_decomposition[['joint_scores']] <- joint_scores
  jive_decomposition[['joint_rank']] <- joint_rank

  jive_decomposition[['joint_rank_sel']] <- joint_rank_sel_results
  class(jive_decomposition) <- "rajive"
  jive_decomposition
}


#' The singular value threshold.
#'
#' Computes the singular value threshold for the data matrix (half way between the rank and rank + 1 singluar value).
#'
#' @param singular_values Numeric. The singular values.
#' @param rank Integer. The rank of the approximation.
#'
#' @details
#' When \code{rank < length(singular_values)} the threshold is the midpoint
#' between \code{singular_values[rank]} and \code{singular_values[rank + 1]},
#' which is the standard gap-midpoint rule.
#'
#' \strong{Boundary heuristic:} when \code{rank == length(singular_values)}
#' (i.e. the requested rank equals the data rank so there is no \code{rank+1}
#' singular value) the function returns \code{0.5 * singular_values[rank]}.
#' This is an arbitrary heuristic in the interval \code{(0, sv[rank])}.  If
#' you encounter this edge case, verify that \code{initial_signal_rank} is
#' strictly less than \code{min(nrow(X), ncol(X))} for each block.

get_sv_threshold <- function(singular_values, rank){
  # When rank == length(singular_values) there is no rank+1 entry; use half
  # the last singular value (midpoint to the implicit floor at 0) to avoid NA.
  if (rank == length(singular_values)) {
    0.5 * singular_values[rank]
  } else {
    0.5 * (singular_values[rank] + singular_values[rank + 1L])
  }
}


#' Computes the joint scores.
#'
#' Estimate the joint rank with the wedin bound, compute the signal scores SVD, double check each joint component.
#'
#' @param blocks List. A list of the data matrices.
#' @param block_svd List. The SVD of the data blocks.
#' @param initial_signal_ranks Numeric vector. Initial signal ranks estimates.
#' @param sv_thresholds Numeric vector. The singular value thresholds from the initial signal rank estimates.
#' @param n_wedin_samples Integer. Number of wedin bound samples to draw for each data matrix.
#' @param n_rand_dir_samples Integer. Number of random direction bound samples to draw.
#' @param joint_rank Integer or NA. User specified joint_rank. If NA will be estimated from data.
#' @param n_perm_samples Integer or NA. Number of permutation samples for the permutation-based
#'   joint rank threshold. See \code{\link{Rajive}} for full description.
#' @param num_cores Integer. Number of cores for parallel resampling.
#' @param identifiability_norm Character. \code{"l2"} uses Euclidean norm for
#'   the post-selection identifiability filter; \code{"l1"} matches original
#'   RaJIVE's one-column \code{norm(score)} behavior.
#' @importFrom stats quantile
#'
#'
#'

get_joint_scores_robustH <- function(blocks, block_svd, initial_signal_ranks, sv_thresholds,
                                     n_wedin_samples=1000, n_rand_dir_samples=1000,
                                     joint_rank=NA,
                                     n_perm_samples=NA,
                                     num_cores=2,
                                     identifiability_norm=c("l2", "l1")){

  identifiability_norm <- .match_identifiability_norm(identifiability_norm)

  if(is.na(n_wedin_samples) & is.na(n_rand_dir_samples) & is.na(n_perm_samples) & is.na(joint_rank)){
    stop('at least one of n_wedin_samples, n_rand_dir_samples, n_perm_samples, or joint_rank must not be NA',
         call.=FALSE)
  }

  K <- length(blocks)
  n_obs <- dim(blocks[[1]])[1]

  # SVD of the signal scores matrix -----------------------------------------
  signal_scores <- list()
  for(k in 1:K){
    signal_scores[[k]] <- block_svd[[k]][['u']][, 1:initial_signal_ranks[k], drop = FALSE]
  }

  M <- do.call(cbind, signal_scores)
  M_svd <- get_svd_robustH(M, rank=min(initial_signal_ranks))


  # estimate joint rank with wedin bound and random direction bound -------------------------------------------------------------

  rank_sel_results  <- list()
  rank_sel_results[['obs_svals']] <- M_svd[['d']]
  rank_sel_results[['identifiability_norm']] <- identifiability_norm

  if(is.na(joint_rank)){

    # maybe comptue wedin bound
    if(!is.na(n_wedin_samples)){

      # NOTE: use explicit loop (not mapply) so each block gets its own
      # signal_rank correctly, avoiding variable-capture from the outer loop.
      block_wedin_samples <- t(do.call(rbind, lapply(seq_along(blocks), function(k)
        get_wedin_bound_samples(blocks[[k]],
                                block_svd[[k]],
                                signal_rank=initial_signal_ranks[k],
                                num_samples=n_wedin_samples,
                                num_cores=num_cores))))




      wedin_samples <-  K - colSums(block_wedin_samples)
      wedin_svsq_threshold <- quantile(wedin_samples, .05)

      rank_sel_results[['wedin']] <- list(block_wedin_samples=block_wedin_samples,
                                          wedin_samples=wedin_samples,
                                          wedin_svsq_threshold=wedin_svsq_threshold)
    } else{
      wedin_svsq_threshold <- NA
    }

    # maybe compute random direction bound or permutation bound
    # permutation replaces random direction when n_perm_samples is set
    perm_svsq_threshold    <- NA
    rand_dir_svsq_threshold <- NA

    if (!is.na(n_perm_samples)) {

      perm_draws <- get_perm_bound_robustH(block_svd = block_svd,
                                           initial_signal_ranks = initial_signal_ranks,
                                           num_samples = n_perm_samples,
                                           num_cores = num_cores)
      perm_svsq_threshold <- quantile(perm_draws, .95)

      rank_sel_results[['perm']] <- list(perm_samples = perm_draws,
                                         perm_svsq_threshold = perm_svsq_threshold)

    }
    if (!is.na(n_rand_dir_samples)) {

      rand_dir_samples <- get_random_direction_bound_robustH(n_obs = n_obs,
                                                             dims = initial_signal_ranks,
                                                             num_samples = n_rand_dir_samples,
                                                             num_cores = num_cores)
      rand_dir_svsq_threshold <- quantile(rand_dir_samples, .95)

      rank_sel_results[['rand_dir']] <- list(rand_dir_samples = rand_dir_samples,
                                             rand_dir_svsq_threshold = rand_dir_svsq_threshold)
    }

    overall_sv_sq_threshold <- max(wedin_svsq_threshold, rand_dir_svsq_threshold,
                                   perm_svsq_threshold, na.rm = TRUE)
    joint_rank_estimate <- sum(M_svd[['d']]^2 > overall_sv_sq_threshold)

    rank_sel_results[['overall_sv_sq_threshold']] <- overall_sv_sq_threshold
    rank_sel_results[['joint_rank_estimate']] <- joint_rank_estimate


  } else { # user provided joint rank
    joint_rank_estimate <- joint_rank
    rank_sel_results[['joint_rank_estimate']] <- joint_rank
  }


  # estimate joint score space ------------------------------------


  joint_scores <- M_svd[['u']][ , seq_len(joint_rank_estimate), drop=FALSE]

  # reconsider joint score space ------------------------------------
  # remove columns of joint_scores that have a
  # trivial projection from one of the data matrices

  to_remove <- c()
  for(k in seq_len(K)){
    for(j in seq_len(joint_rank_estimate)){

      score <- t(blocks[[k]]) %*% joint_scores[ , j]
      sv <- .identifiability_projection_norm(score, identifiability_norm)

      if(sv < sv_thresholds[[k]]){
        message('removing column ', j)
        to_remove <- c(to_remove, j)
        break
      }
    }

  }
  to_keep <- setdiff(seq_len(joint_rank_estimate), to_remove)
  joint_rank <- length(to_keep)
  joint_scores <- joint_scores[ , to_keep, drop=FALSE]

  rank_sel_results[['identif_dropped']] <- if (length(to_remove) > 0L)
    sort(unique(to_remove)) else integer(0L)

  list(joint_scores=joint_scores, rank_sel_results=rank_sel_results)
}


#' Computes the final JIVE decomposition.
#'
#' Computes X = J + I + E for a single data block and the respective SVDs.
#'
#'
#' @param X Matrix. The original data matrix.
#' @param joint_scores Matrix. The basis of the joint space (dimension n x joint_rank).
#' @param sv_threshold Numeric vector. The singular value thresholds from the initial signal rank estimates.
#' @param full Boolean. Do we compute the full J, I matrices or just svd
#'
get_final_decomposition_robustH <- function(X, joint_scores, sv_threshold, full=TRUE){

  jive_decomposition <- list()
  jive_decomposition[['individual']] <- get_individual_decomposition_robustH(X, joint_scores, sv_threshold, full)
  jive_decomposition[['joint']] <- get_joint_decomposition_robustH(X, joint_scores, full)


  if(full){
    jive_decomposition[['noise']] <- X - (jive_decomposition[['joint']][['full']] +
                                            jive_decomposition[['individual']][['full']])
  } else{
    jive_decomposition[['noise']] <- NA
  }

  jive_decomposition
}


#' Computes the individual matrix for a data block.
#'
#' @param X Matrix. The original data matrix.
#' @param joint_scores Matrix. The basis of the joint space (dimension n x joint_rank).
#' @param sv_threshold Numeric vector. The singular value thresholds from the initial signal rank estimates.
#' @param full Boolean. Do we compute the full J, I matrices or just the SVD (set to FALSE to save memory).

get_individual_decomposition_robustH <- function(X, joint_scores, sv_threshold, full=TRUE){

  # When joint_rank == 0 (after identifiability filtering), there is no
  # joint subspace to project out: X_orthog == X.  Avoid forming an empty
  # n x 0 / 0 x n product and the n x n identity.
  if (is.null(joint_scores) || ncol(joint_scores) == 0L) {
    X_orthog <- X
  } else {
    # Avoid building the full n x n projector: I - U U^T applied to X is
    # X - U (U^T X) with two small mat-muls.
    X_orthog <- X - joint_scores %*% (t(joint_scores) %*% X)
  }

  indiv_decomposition <- get_svd_robustH(X_orthog)

  indiv_rank <- sum(indiv_decomposition[['d']] > sv_threshold)

  indiv_decomposition <- truncate_svd(decomposition=indiv_decomposition,
                                      rank=indiv_rank)
  if(full){
    indiv_decomposition[['full']] <- svd_reconstruction(indiv_decomposition)
  } else{
    indiv_decomposition[['full']] <- NA
  }

  indiv_decomposition[['rank']] <- indiv_rank
  indiv_decomposition
}


#' Computes the individual matrix for a data block
#'
#' @param X Matrix. The original data matrix.
#' @param joint_scores Matrix. The basis of the joint space (dimension n x joint_rank).
#' @param full Boolean. Do we compute the full J, I matrices or just the SVD (set to FALSE to save memory).

get_joint_decomposition_robustH <- function(X, joint_scores, full=TRUE){

  joint_rank <- dim(joint_scores)[2]

  # Degenerate case: the identifiability filter removed every joint
  # component (or the user supplied joint_rank = 0).  Return a zero
  # decomposition so downstream `final_decomposition` can compute the
  # noise as X - 0 - I correctly without invoking the C++ M-estimator
  # on a zero matrix (which fails to converge: `solve(): solution not
  # found`).
  if (is.null(joint_rank) || joint_rank == 0L) {
    n <- dim(X)[1]
    p <- dim(X)[2]
    joint_decomposition <- list(
      u    = matrix(0, nrow = n, ncol = 0L),
      d    = numeric(0L),
      v    = matrix(0, nrow = p, ncol = 0L),
      full = if (isTRUE(full)) matrix(0, nrow = n, ncol = p) else NA,
      rank = 0L
    )
    return(joint_decomposition)
  }

  J <-  joint_scores %*% t(joint_scores) %*% X

  # Exact low-rank decomposition: if J = U U^T X with U orthonormal,
  # then SVD(U^T X) = U_small D V^T implies J = (U U_small) D V^T.
  # This avoids an unnecessary robust M-estimation pass on rank-joint_rank J.
  B <- t(joint_scores) %*% X
  B_svd <- svd(B, nu = joint_rank, nv = joint_rank)
  joint_decomposition <- list(
    u = joint_scores %*% B_svd$u,
    d = B_svd$d,
    v = B_svd$v
  )

  if(full){
    joint_decomposition[['full']] <- J
  } else{
    joint_decomposition[['full']] <- NA
  }

  joint_decomposition[['rank']] <- joint_rank
  joint_decomposition

}
