# Internal bootstrap engine for RaJIVE inference helpers.

.bootstrap_resample_indices <- function(n, sample_frac = 0.8,
                                        cluster = NULL,
                                        strata = NULL,
                                        resample = NULL) {
  if (!is.numeric(n) || length(n) != 1L || n < 2L) {
    cli::cli_abort("`n` must be a sample count of at least 2.")
  }
  n <- as.integer(n)
  if (!is.numeric(sample_frac) || length(sample_frac) != 1L ||
      !is.finite(sample_frac) || sample_frac <= 0 || sample_frac > 1) {
    cli::cli_abort("`sample_frac` must be in (0, 1].")
  }

  if (is.null(resample)) {
    resample <- if (is.null(cluster)) "observation" else "cluster"
  }
  resample <- match.arg(resample, c("observation", "cluster"))

  if (resample == "observation") {
    size <- max(2L, floor(sample_frac * n))
    return(sample.int(n, size = size, replace = TRUE))
  }

  if (is.null(cluster)) {
    cli::cli_abort("`cluster` must be supplied when `resample = 'cluster'`.")
  }
  if (length(cluster) != n) {
    cli::cli_abort("`cluster` must have one value per sample.")
  }
  if (!is.null(strata) && length(strata) != n) {
    cli::cli_abort("`strata` must have one value per sample.")
  }

  cluster_f <- factor(cluster)
  row_by_cluster <- split(seq_len(n), cluster_f)

  if (is.null(strata)) {
    cluster_groups <- list(all = names(row_by_cluster))
  } else {
    strata_f <- factor(strata)
    cluster_strata <- vapply(row_by_cluster, function(idx) {
      vals <- unique(as.character(strata_f[idx]))
      if (length(vals) != 1L) {
        cli::cli_abort("Each cluster must belong to exactly one stratum.")
      }
      vals
    }, character(1L))
    cluster_groups <- split(names(cluster_strata), cluster_strata)
  }

  sampled <- unlist(lapply(cluster_groups, function(ids) {
    size <- max(1L, floor(sample_frac * length(ids)))
    draw <- sample(ids, size = size, replace = TRUE)
    unlist(row_by_cluster[draw], use.names = FALSE)
  }), use.names = FALSE)

  as.integer(sampled)
}

.scatter_scores_to_reference_rows <- function(scores, idx, n_ref) {
  out <- matrix(NA_real_, nrow = n_ref, ncol = ncol(scores))
  for (i in unique(idx)) {
    pos <- which(idx == i)
    out[i, ] <- colMeans(scores[pos, , drop = FALSE])
  }
  out
}

.validate_bootstrap_score_request <- function(score_type, score_block, K) {
  score_type <- match.arg(score_type, c("joint", "individual"))
  if (score_type == "individual") {
    if (is.null(score_block)) {
      cli::cli_abort("`score_block` must be supplied for individual score bootstrap.")
    }
    if (!is.numeric(score_block) || length(score_block) != 1L ||
        is.na(score_block) || score_block < 1L || score_block > K) {
      cli::cli_abort("`score_block` must be an integer block index between 1 and {K}.")
    }
  }
  invisible(TRUE)
}

.extract_score_matrix <- function(fit, score_type = c("joint", "individual"),
                                  score_block = NULL) {
  score_type <- match.arg(score_type)
  if (is.null(fit)) return(NULL)
  if (score_type == "joint") {
    return(fit$joint_scores)
  }

  if (is.null(score_block) || is.null(fit$block_decomps)) return(NULL)
  idx <- 3L * (as.integer(score_block) - 1L) + 1L
  if (idx < 1L || idx > length(fit$block_decomps)) return(NULL)
  fit$block_decomps[[idx]]$u
}

.align_and_scatter_scores <- function(ref_scores, bs_scores, idx, n_ref,
                                      align_to = c("reference", "none")) {
  align_to <- match.arg(align_to)
  if (is.null(ref_scores) || is.null(bs_scores)) return(NULL)
  if (is.vector(ref_scores)) ref_scores <- matrix(ref_scores, ncol = 1L)
  if (is.vector(bs_scores)) bs_scores <- matrix(bs_scores, ncol = 1L)
  if (ncol(ref_scores) == 0L || ncol(bs_scores) == 0L) return(NULL)

  n_use <- min(ncol(ref_scores), ncol(bs_scores))
  if (n_use < 1L || nrow(bs_scores) != length(idx)) return(NULL)

  ref_sub <- ref_scores[idx, seq_len(n_use), drop = FALSE]
  bs_sub <- bs_scores[, seq_len(n_use), drop = FALSE]
  ok <- stats::complete.cases(ref_sub, bs_sub)
  if (sum(ok) < max(3L, n_use)) return(NULL)

  if (align_to == "reference") {
    bs_aligned <- tryCatch({
      Q <- .procrustes_align(ref_sub[ok, , drop = FALSE],
                             bs_sub[ok, , drop = FALSE])
      bs_sub %*% Q
    }, error = function(e) NULL)
    if (is.null(bs_aligned) || any(!is.finite(bs_aligned[ok, , drop = FALSE]))) {
      return(NULL)
    }
  } else {
    bs_aligned <- bs_sub
  }

  list(
    scores = .scatter_scores_to_reference_rows(bs_aligned, idx, n_ref),
    sample_scores = bs_aligned,
    n_use = n_use
  )
}

.bootstrap_score_targets <- function(ajive_output, targets, B, n_ref,
                                     require_positive = TRUE) {
  out <- vector("list", nrow(targets))
  names(out) <- targets$key
  for (i in seq_len(nrow(targets))) {
    score_block <- if (is.na(targets$block[[i]])) NULL else targets$block[[i]]
    ref <- .extract_score_matrix(ajive_output, targets$source[[i]], score_block)
    if (is.vector(ref)) ref <- matrix(ref, ncol = 1L)
    if (is.null(ref) || ncol(ref) == 0L) {
      if (require_positive) {
        cli::cli_abort("Requested {.val {targets$source[[i]]}} score target has zero fitted rank.")
      }
      next
    }
    out[[i]] <- list(
      source = targets$source[[i]],
      block = targets$block[[i]],
      ref_scores = ref,
      scores = array(NA_real_, dim = c(n_ref, ncol(ref), B))
    )
  }
  out
}

.component_var_explained <- function(fit, blocks, n_comp) {
  K <- length(blocks)
  out <- matrix(NA_real_, nrow = K, ncol = n_comp)
  for (k in seq_len(K)) {
    decomp <- fit$block_decomps[[3L * (k - 1L) + 2L]]
    d <- decomp[["d"]]
    if (length(d) == 0L) next
    n_use <- min(length(d), n_comp)
    denom <- norm(blocks[[k]], type = "F")^2
    out[k, seq_len(n_use)] <- d[seq_len(n_use)]^2 / denom
  }
  out
}

.rajive_bootstrap <- function(ajive_output, blocks, initial_signal_ranks,
                              B = 100L,
                              sample_frac = 0.8,
                              cluster = NULL,
                              strata = NULL,
                              resample = NULL,
                              align_to = c("reference", "first_replicate", "none"),
                              num_cores = 1L,
                              keep = c("loadings", "scores", "joint_rank",
                                       "component_cors", "indices",
                                       "var_explained"),
                              score_type = c("joint", "individual"),
                              score_block = NULL,
                              ...) {
  align_to <- match.arg(align_to)
  score_type <- match.arg(score_type)
  keep <- match.arg(keep, c("loadings", "scores", "joint_rank",
                            "component_cors", "indices", "var_explained"),
                    several.ok = TRUE)

  .validate_feature_space(blocks)
  .validate_matched_samples(blocks)
  if (length(initial_signal_ranks) != length(blocks)) {
    cli::cli_abort("`initial_signal_ranks` must have one entry per block.")
  }
  B <- as.integer(B)
  if (B < 1L) {
    cli::cli_abort("`B` must be a positive integer.")
  }

  K <- length(blocks)
  n_ref <- nrow(blocks[[1L]])
  .validate_bootstrap_score_request(score_type, score_block, K)
  if (score_type != "joint" && any(keep %in% c("loadings", "var_explained"))) {
    cli::cli_abort("Individual score bootstrap currently supports score payloads only.")
  }
  needs_reference <- any(keep %in% c("loadings", "scores", "component_cors",
                                     "var_explained"))
  n_comp <- 0L
  ref_scores <- NULL
  if (needs_reference) {
    if (is.null(ajive_output) || !inherits(ajive_output, "rajive")) {
      cli::cli_abort("`ajive_output` must be a fitted rajive object for the requested bootstrap payload.")
    }
    ref_scores <- .extract_score_matrix(ajive_output, score_type, score_block)
    n_comp <- if (!is.null(ref_scores)) ncol(ref_scores) else 0L
    if (is.null(n_comp) || is.na(n_comp) || n_comp < 1L) {
      cli::cli_abort("`ajive_output` must contain at least one requested score component.")
    }
  }

  ref_loadings <- NULL
  if ("loadings" %in% keep && !is.null(ajive_output$block_decomps)) {
    ref_loadings <- lapply(seq_len(K), function(k) {
      ajive_output$block_decomps[[3L * (k - 1L) + 2L]]$v
    })
  }

  out <- list()
  if ("loadings" %in% keep) {
    out$loadings <- lapply(seq_len(K), function(k) {
      array(NA_real_, dim = c(ncol(blocks[[k]]), n_comp, B))
    })
    names(out$loadings) <- .default_block_names(blocks)
  }
  if ("scores" %in% keep) {
    out$scores <- array(NA_real_, dim = c(n_ref, n_comp, B))
  }
  if ("joint_rank" %in% keep) out$joint_rank <- rep(NA_integer_, B)
  if ("component_cors" %in% keep) {
    out$component_cors <- matrix(NA_real_, nrow = B, ncol = n_comp)
  }
  if ("indices" %in% keep) out$indices <- vector("list", B)
  if ("var_explained" %in% keep) {
    out$var_explained <- array(NA_real_, dim = c(K, n_comp, B))
  }

  dots <- list(...)
  for (b in seq_len(B)) {
    idx <- .bootstrap_resample_indices(
      n = n_ref,
      sample_frac = sample_frac,
      cluster = cluster,
      strata = strata,
      resample = resample
    )
    b_list <- lapply(blocks, function(x) x[idx, , drop = FALSE])
    fit_b <- tryCatch(
      do.call(Rajive, c(list(b_list, initial_signal_ranks),
                        dots)),
      error = function(e) NULL
    )

    if ("indices" %in% keep) out$indices[[b]] <- idx
    if (is.null(fit_b)) next
    if ("joint_rank" %in% keep) out$joint_rank[[b]] <- fit_b$joint_rank

    bs <- .extract_score_matrix(fit_b, score_type, score_block)
    if (any(keep %in% c("scores", "component_cors")) &&
        !is.null(bs) && ncol(bs) > 0L && !is.null(ref_scores)) {
      aligned <- .align_and_scatter_scores(
        ref_scores = ref_scores,
        bs_scores = bs,
        idx = idx,
        n_ref = n_ref,
        align_to = if (align_to == "reference") "reference" else "none"
      )
      if (!is.null(aligned)) {
        n_use <- aligned$n_use
        if ("scores" %in% keep) {
          out$scores[, seq_len(n_use), b] <- aligned$scores
        }
        if ("component_cors" %in% keep) {
          for (j in seq_len(n_use)) {
            out$component_cors[b, j] <- suppressWarnings(abs(stats::cor(
              ref_scores[idx, j], aligned$sample_scores[, j],
              use = "pairwise.complete.obs"
            )))
          }
        }
      }
    }

    if ("loadings" %in% keep && !is.null(fit_b$block_decomps) &&
        !is.null(ref_loadings)) {
      for (k in seq_len(K)) {
        L_b <- fit_b$block_decomps[[3L * (k - 1L) + 2L]]$v
        if (is.null(L_b) || ncol(L_b) == 0L) next
        n_use <- min(ncol(L_b), n_comp)
        L_sub <- L_b[, seq_len(n_use), drop = FALSE]
        if (align_to == "reference") {
          Q <- .procrustes_align(ref_loadings[[k]][, seq_len(n_use), drop = FALSE], L_sub)
          L_sub <- L_sub %*% Q
        }
        out$loadings[[k]][, seq_len(n_use), b] <- L_sub
      }
    }

    if ("var_explained" %in% keep && !is.null(fit_b$block_decomps)) {
      out$var_explained[, , b] <- .component_var_explained(fit_b, b_list, n_comp)
    }
  }

  out[keep]
}

.rajive_bootstrap_scores_multi <- function(ajive_output, blocks,
                                           initial_signal_ranks, targets,
                                           B = 100L, sample_frac = 0.8,
                                           cluster = NULL, strata = NULL,
                                           resample = NULL,
                                           align_to = c("reference", "none"),
                                           num_cores = 1L, ...) {
  align_to <- match.arg(align_to)
  .validate_feature_space(blocks)
  .validate_matched_samples(blocks)
  if (!inherits(ajive_output, "rajive")) {
    cli::cli_abort("`ajive_output` must be a fitted rajive object.")
  }
  if (length(initial_signal_ranks) != length(blocks)) {
    cli::cli_abort("`initial_signal_ranks` must have one entry per block.")
  }
  B <- as.integer(B)
  if (B < 1L) cli::cli_abort("`B` must be a positive integer.")

  K <- length(blocks)
  n_ref <- nrow(blocks[[1L]])
  for (i in seq_len(nrow(targets))) {
    block_i <- if (is.na(targets$block[[i]])) NULL else targets$block[[i]]
    .validate_bootstrap_score_request(targets$source[[i]], block_i, K)
  }
  out <- .bootstrap_score_targets(ajive_output, targets, B, n_ref)

  dots <- list(...)
  for (b in seq_len(B)) {
    idx <- .bootstrap_resample_indices(
      n = n_ref,
      sample_frac = sample_frac,
      cluster = cluster,
      strata = strata,
      resample = resample
    )
    b_list <- lapply(blocks, function(x) x[idx, , drop = FALSE])
    fit_b <- tryCatch(
      do.call(Rajive, c(list(b_list, initial_signal_ranks), dots)),
      error = function(e) NULL
    )
    if (is.null(fit_b)) next

    for (i in seq_along(out)) {
      target <- out[[i]]
      if (is.null(target)) next
      score_block <- if (is.na(target$block)) NULL else target$block
      bs <- .extract_score_matrix(fit_b, target$source, score_block)
      aligned <- .align_and_scatter_scores(
        ref_scores = target$ref_scores,
        bs_scores = bs,
        idx = idx,
        n_ref = n_ref,
        align_to = align_to
      )
      if (is.null(aligned)) next
      out[[i]]$scores[, seq_len(aligned$n_use), b] <- aligned$scores
    }
  }

  lapply(out, function(x) list(scores = x$scores))
}
