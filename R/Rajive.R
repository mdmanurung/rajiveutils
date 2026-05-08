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


#' Robust Angle-based Joint and Individual Variation Explained (RaJIVE)
#'
#' Computes the robust aJIVE decomposition of a list of multi-view data
#' matrices into joint, block-individual, and residual (noise) components.
#' The robust SVD step uses an M-estimator (Huber loss) so the decomposition
#' is resistant to a moderate fraction of element-wise outliers in any block.
#'
#' Joint rank is selected by comparing the squared common singular values of
#' the stacked signal-score matrix to a threshold derived from a Wedin bound
#' and either a random-direction bound (default) or a permutation bound
#' (when \code{n_perm_samples} is supplied). The combined threshold is
#' \code{max(wedin, rand_dir | perm)}; this heuristic is conservative in
#' practice but does not carry a formal FWER/FDR guarantee for rank
#' selection (see \code{StatisticalAudits.md}, Finding 6).
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
#'   bound samples to draw. Default \code{1000}. Ignored when
#'   \code{n_perm_samples} is supplied.
#' @param joint_rank Integer or \code{NA}. User-specified joint rank. When
#'   \code{NA} (default) the joint rank is estimated from the data.
#' @param n_perm_samples Integer or \code{NA}. Number of permutation samples
#'   for the permutation-based joint-rank threshold. When supplied, replaces
#'   the random-direction bound: rows of each block's signal scores are
#'   independently shuffled, the leading squared singular value of the
#'   stacked permuted matrix is recorded, and its 95th percentile is used
#'   as the threshold. Permutation takes precedence when both are set.
#'   Default \code{NA}.
#' @param num_cores Positive integer. Number of cores to use for parallel
#'   block SVDs, Wedin / random-direction / permutation bound resampling.
#'   Default \code{1L} (serial). When \code{> 1} the function transparently
#'   uses \code{\link[parallel]{mclapply}} on Unix and a
#'   \code{\link[parallel]{makePSOCKcluster}} on Windows. The RNG kind is
#'   temporarily set to \code{"L'Ecuyer-CMRG"} for reproducible parallel
#'   sampling (with a one-time warning); set
#'   \code{RNGkind("L'Ecuyer-CMRG")} explicitly to silence it.
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
#'       components dropped by the L2 identifiability filter.}
#'   }
#'
#' @section Reproducibility:
#' For exactly reproducible results, set \code{set.seed()} (and, for
#' \code{num_cores > 1}, \code{RNGkind("L'Ecuyer-CMRG")}) before calling
#' \code{Rajive()}.
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
                           num_cores=1L)
{

  num_cores <- max(1L, as.integer(num_cores))

  # Phase 4 boundary validation: abort early with clear classes/messages.
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

  degenerate_blocks <- which(vapply(blocks, function(x)
    any(apply(x, 2L, stats::sd) < .Machine$double.eps^0.5), logical(1L)))
  if (length(degenerate_blocks) > 0L) {
    cli::cli_abort(
      c("Degenerate block detected: one or more columns have near-zero variance.",
        "x" = "Block index(es): {.val {degenerate_blocks}}"),
      class = "rajiveplus_degenerate_block"
    )
  }

  n_obs <- nrow(blocks[[1L]])
  if (n_obs < sum(initial_signal_ranks)) {
    cli::cli_warn(
      c("`n` is smaller than sum(initial_signal_ranks); estimation may be underdetermined.",
        "i" = "n = {.val {n_obs}}, sum(initial_signal_ranks) = {.val {sum(initial_signal_ranks)}}"),
      class = "rajiveplus_underdetermined"
    )
  }

  # Phase 2 reproducibility: parallel paths require L'Ecuyer-CMRG.
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
                                  num_cores=num_cores)
  joint_rank_sel_results <- out$rank_sel_results
  joint_scores <- out$joint_scores

  joint_rank <- dim(joint_scores)[2]


  # step 3: final decomposition -----------------------------------------------------

  block_decomps <- mapply(function(l,m)
    get_final_decomposition_robustH(l, joint_scores = joint_scores, m), blocks,
    sv_thresholds )



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

get_sv_threshold <- function(singular_values, rank){
  # W-H1: boundary guard — when rank == length(singular_values) there is no
  # rank+1 entry; use half the last singular value (midpoint to the implicit
  # floor at 0) to avoid returning NA.
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
#' @importFrom stats quantile
#'
#'
#'

get_joint_scores_robustH <- function(blocks, block_svd, initial_signal_ranks, sv_thresholds,
                                     n_wedin_samples=1000, n_rand_dir_samples=1000,
                                     joint_rank=NA,
                                     n_perm_samples=NA,
                                     num_cores=2){


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

    } else if (!is.na(n_rand_dir_samples)) {

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


  joint_scores <- M_svd[['u']][ , 1:joint_rank_estimate, drop=FALSE]

  # reconsider joint score space ------------------------------------
  # remove columns of joint_scores that have a
  # trivial projection from one of the data matrices

  to_remove <- c()
  for(k in 1:K){
    for(j in 1:joint_rank_estimate){

      score <- t(blocks[[k]]) %*% joint_scores[ , j]
      # Spectral / L2 norm to match the scale of `sv_thresholds`, which are
      # derived from singular values.  `norm()` defaults to type = "O"
      # (max absolute column sum = sum(|score|) for a one-column matrix),
      # which is systematically larger than the L2 norm and biased the
      # identifiability filter toward keeping components.  See audit
      # finding #3.
      sv <- sqrt(sum(score^2))

      if(sv < sv_thresholds[[k]]){
        message('removing column ', j)
        to_remove <- c(to_remove, j)
        break
      }
    }

  }
  to_keep <- setdiff(1:joint_rank_estimate, to_remove)
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

  # Audit perf #14: J = U U^T X is exactly rank `joint_rank` by construction
  # (joint_scores is an orthonormal basis), so plain base::svd() is
  # mathematically equivalent to the M-estimator robust SVD here but ~10-100x
  # faster.  No outliers can survive the projection onto the joint subspace.
  joint_decomposition <- svd(J, nu = joint_rank, nv = joint_rank)
  joint_decomposition <- list(
    u = joint_decomposition$u[, seq_len(joint_rank), drop = FALSE],
    d = joint_decomposition$d[seq_len(joint_rank)],
    v = joint_decomposition$v[, seq_len(joint_rank), drop = FALSE]
  )

  if(full){
    joint_decomposition[['full']] <- J
  } else{
    joint_decomposition[['full']] <- NA
  }

  joint_decomposition[['rank']] <- joint_rank
  joint_decomposition

}



