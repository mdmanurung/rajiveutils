# Native missing-data helpers -------------------------------------------------

.missing_block_names <- function(blocks) {
  if (!is.null(names(blocks)) && all(nzchar(names(blocks)))) {
    names(blocks)
  } else {
    paste0("block", seq_along(blocks))
  }
}

.validate_missing_blocks <- function(blocks) {
  if (!is.list(blocks) || length(blocks) < 1L) {
    cli::cli_abort(
      "`blocks` must be a non-empty list of numeric matrices.",
      class = "rajiveplus_invalid_input"
    )
  }
  if (!all(vapply(blocks, is.matrix, logical(1L)))) {
    cli::cli_abort(
      "Every element of `blocks` must be a matrix.",
      class = "rajiveplus_invalid_input"
    )
  }
  if (!all(vapply(blocks, is.numeric, logical(1L)))) {
    cli::cli_abort(
      "Every element of `blocks` must be numeric.",
      class = "rajiveplus_invalid_input"
    )
  }
  bad_dim <- vapply(blocks, function(x) any(dim(x) == 0L), logical(1L))
  if (any(bad_dim)) {
    cli::cli_abort(
      "Every block must have at least one sample and one feature.",
      class = "rajiveplus_invalid_input"
    )
  }
  n_obs <- vapply(blocks, nrow, integer(1L))
  if (length(unique(n_obs)) != 1L) {
    cli::cli_abort(
      "All blocks must have the same number of samples.",
      class = "rajiveplus_invalid_input"
    )
  }
  invisible(blocks)
}

.normalize_missing_mask <- function(blocks, mask = NULL) {
  .validate_missing_blocks(blocks)
  block_names <- .missing_block_names(blocks)
  names(blocks) <- block_names

  if (is.null(mask)) {
    mask <- lapply(blocks, is.finite)
  } else {
    if (is.matrix(mask) && length(blocks) == 1L) {
      mask <- list(mask)
    }
    if (!is.list(mask) || length(mask) != length(blocks)) {
      cli::cli_abort(
        "`mask` must be a list with one logical matrix per block.",
        class = "rajiveplus_invalid_input"
      )
    }
    names(mask) <- block_names
  }

  for (k in seq_along(blocks)) {
    if (!is.matrix(mask[[k]]) || !identical(dim(mask[[k]]), dim(blocks[[k]]))) {
      cli::cli_abort(
        c("Each `mask` entry must have the same dimensions as its block.",
          "i" = "{block_names[[k]]}: block is {paste(dim(blocks[[k]]), collapse = ' x ')}, mask is {paste(dim(mask[[k]]), collapse = ' x ')}."),
        class = "rajiveplus_invalid_input"
      )
    }
    if (anyNA(mask[[k]])) {
      cli::cli_abort(
        "`mask` entries must not contain missing values.",
        class = "rajiveplus_invalid_input"
      )
    }
    mask[[k]] <- matrix(as.logical(mask[[k]]),
                        nrow = nrow(mask[[k]]),
                        ncol = ncol(mask[[k]]),
                        dimnames = dimnames(blocks[[k]]))
    nonfinite_observed <- mask[[k]] & !is.finite(blocks[[k]])
    if (any(nonfinite_observed)) {
      first <- which(nonfinite_observed, arr.ind = TRUE)[1L, , drop = TRUE]
      cli::cli_abort(
        c("Observed entries must be finite.",
          "i" = "{block_names[[k]]}[{first[[1L]]}, {first[[2L]]}] is non-finite but `mask` marks it observed."),
        class = "rajiveplus_invalid_input"
      )
    }
  }

  patterns <- .classify_missingness(blocks, mask)
  structure(
    list(
      blocks = blocks,
      mask = mask,
      observed = mask,
      patterns = patterns,
      block_names = block_names
    ),
    class = "rajive_missing_mask"
  )
}

.classify_missingness <- function(blocks, mask) {
  .validate_missing_blocks(blocks)
  block_names <- .missing_block_names(blocks)
  names(blocks) <- block_names
  names(mask) <- block_names

  row_ids <- rownames(blocks[[1L]])
  if (is.null(row_ids)) row_ids <- as.character(seq_len(nrow(blocks[[1L]])))

  sample_block_rows <- do.call(rbind, lapply(seq_along(blocks), function(k) {
    obs <- mask[[k]]
    observed_cells <- rowSums(obs)
    data.frame(
      block = block_names[[k]],
      sample = seq_len(nrow(obs)),
      sample_id = row_ids,
      observed_cells = as.integer(observed_cells),
      missing_cells = as.integer(ncol(obs) - observed_cells),
      whole_row_missing = observed_cells == 0L,
      partial_row_missing = observed_cells > 0L & observed_cells < ncol(obs),
      complete_row = observed_cells == ncol(obs),
      stringsAsFactors = FALSE
    )
  }))
  rownames(sample_block_rows) <- NULL

  features <- do.call(rbind, lapply(seq_along(blocks), function(k) {
    obs <- mask[[k]]
    col_ids <- colnames(blocks[[k]])
    if (is.null(col_ids)) col_ids <- as.character(seq_len(ncol(obs)))
    observed_cells <- colSums(obs)
    data.frame(
      block = block_names[[k]],
      feature = seq_len(ncol(obs)),
      feature_id = col_ids,
      observed_cells = as.integer(observed_cells),
      missing_cells = as.integer(nrow(obs) - observed_cells),
      all_missing = observed_cells == 0L,
      partial_missing = observed_cells > 0L & observed_cells < nrow(obs),
      complete_feature = observed_cells == nrow(obs),
      stringsAsFactors = FALSE
    )
  }))
  rownames(features) <- NULL

  cells <- do.call(rbind, lapply(seq_along(blocks), function(k) {
    obs <- mask[[k]]
    idx <- expand.grid(row = seq_len(nrow(obs)),
                       col = seq_len(ncol(obs)),
                       KEEP.OUT.ATTRS = FALSE)
    whole_row <- rowSums(obs) == 0L
    all_feature <- colSums(obs) == 0L
    observed <- obs[cbind(idx$row, idx$col)]
    data.frame(
      block = block_names[[k]],
      row = idx$row,
      col = idx$col,
      observed = observed,
      missing = !observed,
      whole_row_missing = whole_row[idx$row],
      all_feature_missing = all_feature[idx$col],
      scattered_missing = !observed & !whole_row[idx$row] & !all_feature[idx$col],
      stringsAsFactors = FALSE
    )
  }))
  rownames(cells) <- NULL

  block_summary <- do.call(rbind, lapply(seq_along(blocks), function(k) {
    obs <- mask[[k]]
    data.frame(
      block = block_names[[k]],
      nrow = nrow(obs),
      ncol = ncol(obs),
      observed_cells = as.integer(sum(obs)),
      missing_cells = as.integer(sum(!obs)),
      observed_fraction = mean(obs),
      whole_sample_block_rows = sum(rowSums(obs) == 0L),
      all_missing_features = sum(colSums(obs) == 0L),
      stringsAsFactors = FALSE
    )
  }))
  rownames(block_summary) <- NULL

  sample_observed <- Reduce(`+`, lapply(mask, rowSums))
  sample_summary <- data.frame(
    sample = seq_along(sample_observed),
    sample_id = row_ids,
    observed_cells = as.integer(sample_observed),
    missing_cells = as.integer(sum(vapply(mask, ncol, integer(1L))) - sample_observed),
    all_blocks_missing = sample_observed == 0L,
    stringsAsFactors = FALSE
  )

  list(
    sample_block_rows = sample_block_rows,
    features = features,
    cells = cells,
    block_summary = block_summary,
    sample_summary = sample_summary,
    counts = c(
      missing_cells = sum(!unlist(mask, use.names = FALSE)),
      whole_sample_block_rows = sum(sample_block_rows$whole_row_missing),
      scattered_cells = sum(cells$scattered_missing),
      all_missing_features = sum(features$all_missing),
      samples_all_blocks_missing = sum(sample_summary$all_blocks_missing)
    )
  )
}

.validate_native_missing_inputs <- function(blocks, mask, control = NULL) {
  normalized <- .normalize_missing_mask(blocks, mask)
  block_names <- normalized$block_names
  mask <- normalized$mask
  patterns <- normalized$patterns

  empty_blocks <- which(vapply(mask, function(x) !any(x), logical(1L)))
  if (length(empty_blocks) > 0L) {
    cli::cli_abort(
      c("Native missing-data mode cannot fit a fully missing block.",
        "x" = "Fully missing block(s): {paste(block_names[empty_blocks], collapse = ', ')}."),
      class = "rajiveplus_all_missing_block"
    )
  }

  empty_samples <- patterns$sample_summary$sample_id[
    patterns$sample_summary$all_blocks_missing
  ]
  if (length(empty_samples) > 0L) {
    cli::cli_abort(
      c("Native missing-data mode cannot fit samples missing from every block.",
        "x" = "All-missing sample(s): {paste(empty_samples, collapse = ', ')}."),
      class = "rajiveplus_sample_all_missing"
    )
  }

  empty_features <- patterns$features[patterns$features$all_missing, , drop = FALSE]
  if (nrow(empty_features) > 0L) {
    first <- empty_features[1L, , drop = FALSE]
    cli::cli_abort(
      c("Native missing-data mode cannot fit fully missing features.",
        "x" = "First fully missing feature: {first$block} feature {first$feature_id}."),
      class = "rajiveplus_feature_all_missing"
    )
  }

  normalized
}

.validate_native_initial_signal_ranks <- function(initial_signal_ranks,
                                                  n_blocks) {
  if (length(initial_signal_ranks) != n_blocks) {
    cli::cli_abort(
      c("`initial_signal_ranks` must have one entry per block.",
        "x" = "Got {length(initial_signal_ranks)} rank value(s) for {n_blocks} block(s)."),
      class = "rajiveplus_invalid_input"
    )
  }
  invisible(initial_signal_ranks)
}

.center_scale_observed <- function(blocks, mask, center = TRUE, scale = FALSE,
                                   normalize = TRUE) {
  normalized <- .normalize_missing_mask(blocks, mask)
  blocks <- normalized$blocks
  mask <- normalized$mask
  block_names <- normalized$block_names

  out_blocks <- vector("list", length(blocks))
  centers <- scales <- observed_counts <- vector("list", length(blocks))
  frob_norms <- vector("list", length(blocks))
  names(out_blocks) <- names(centers) <- names(scales) <-
    names(observed_counts) <- names(frob_norms) <- block_names

  for (k in seq_along(blocks)) {
    x <- blocks[[k]]
    obs <- mask[[k]]
    n_obs <- colSums(obs)
    observed_counts[[k]] <- as.integer(n_obs)

    center_vec <- numeric(ncol(x))
    if (isTRUE(center)) {
      center_vec <- vapply(seq_len(ncol(x)), function(j) {
        mean(x[obs[, j], j], na.rm = TRUE)
      }, numeric(1L))
    }
    center_vec[!is.finite(center_vec)] <- 0

    scale_vec <- rep(1, ncol(x))
    if (isTRUE(scale)) {
      scale_vec <- vapply(seq_len(ncol(x)), function(j) {
        vals <- x[obs[, j], j] - center_vec[[j]]
        if (length(vals) <= 1L) {
          1
        } else {
          sqrt(sum(vals^2) / (length(vals) - 1L))
        }
      }, numeric(1L))
      scale_vec[!is.finite(scale_vec) | scale_vec <= .Machine$double.eps] <- 1
    }

    z <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    transformed <- sweep(sweep(x, 2L, center_vec, "-"), 2L, scale_vec, "/")
    z[obs] <- transformed[obs]

    frob <- 1
    if (isTRUE(normalize)) {
      frob <- sqrt(sum(z[obs]^2))
      if (!is.finite(frob) || frob <= .Machine$double.eps) {
        frob <- 1
      }
      z[obs] <- z[obs] / frob
    }

    out_blocks[[k]] <- z
    centers[[k]] <- center_vec
    scales[[k]] <- scale_vec
    frob_norms[[k]] <- frob
  }

  list(
    blocks = out_blocks,
    mask = mask,
    center = centers,
    scale = scales,
    frob_norm = frob_norms,
    observed_counts = observed_counts,
    center_enabled = isTRUE(center),
    scale_enabled = isTRUE(scale),
    normalize_enabled = isTRUE(normalize),
    patterns = normalized$patterns,
    block_names = block_names
  )
}

.backtransform_reconstruction <- function(recon, preprocess) {
  if (!is.list(recon) || length(recon) != length(preprocess$center)) {
    cli::cli_abort(
      "`recon` must be a list with one matrix per preprocessed block.",
      class = "rajiveplus_invalid_input"
    )
  }

  block_names <- names(preprocess$center)
  names(recon) <- block_names
  out <- vector("list", length(recon))
  names(out) <- block_names

  for (k in seq_along(recon)) {
    z <- recon[[k]]
    if (!is.matrix(z)) {
      cli::cli_abort(
        "Every reconstruction must be a matrix.",
        class = "rajiveplus_invalid_input"
      )
    }
    if (!identical(dim(z), dim(preprocess$mask[[k]]))) {
      cli::cli_abort(
        "Reconstruction dimensions must match the preprocessing mask.",
        class = "rajiveplus_invalid_input"
      )
    }
    x <- z * preprocess$frob_norm[[k]]
    x <- sweep(x, 2L, preprocess$scale[[k]], "*")
    x <- sweep(x, 2L, preprocess$center[[k]], "+")
    dimnames(x) <- dimnames(z)
    out[[k]] <- x
  }

  out
}

.masked_frob_norm <- function(x, mask) {
  if (!identical(dim(x), dim(mask))) {
    cli::cli_abort(
      "`x` and `mask` must have the same dimensions.",
      class = "rajiveplus_invalid_input"
    )
  }
  sqrt(sum(x[mask]^2))
}

.masked_crossprod <- function(x, y, mask) {
  if (!identical(dim(x), dim(y)) || !identical(dim(x), dim(mask))) {
    cli::cli_abort(
      "`x`, `y`, and `mask` must have the same dimensions.",
      class = "rajiveplus_invalid_input"
    )
  }
  sum(x[mask] * y[mask])
}

.masked_residual_ss <- function(x, fitted, mask) {
  if (!identical(dim(x), dim(fitted)) || !identical(dim(x), dim(mask))) {
    cli::cli_abort(
      "`x`, `fitted`, and `mask` must have the same dimensions.",
      class = "rajiveplus_invalid_input"
    )
  }
  sum((x[mask] - fitted[mask])^2)
}

.masked_variance_explained <- function(x, joint, individual, mask) {
  if (!identical(dim(x), dim(joint)) ||
      !identical(dim(x), dim(individual)) ||
      !identical(dim(x), dim(mask))) {
    cli::cli_abort(
      "`x`, `joint`, `individual`, and `mask` must have the same dimensions.",
      class = "rajiveplus_invalid_input"
    )
  }
  denom <- sum(x[mask]^2)
  if (!is.finite(denom) || denom <= .Machine$double.eps) denom <- 1
  joint_prop <- sum(joint[mask]^2) / denom
  individual_prop <- sum(individual[mask]^2) / denom
  residual_prop <- 1 - joint_prop - individual_prop
  c(
    Joint = joint_prop,
    Indiv = individual_prop,
    Resid = residual_prop,
    observed_fraction = mean(mask),
    observed_ss = denom
  )
}

.masked_orthogonality <- function(a, b, mask) {
  inner <- .masked_crossprod(a, b, mask)
  denom <- .masked_frob_norm(a, mask) * .masked_frob_norm(b, mask)
  scaled <- if (denom <= .Machine$double.eps) 0 else abs(inner / denom)
  list(
    crossprod = inner,
    max_abs_crossprod = scaled,
    observed_fraction = mean(mask)
  )
}

#' Native missing-data control options
#'
#' Constructs a control list for \code{\link{Rajive}} when
#' \code{missing = "native"}.
#'
#' @param center,scale,normalize Logical preprocessing flags applied with
#'   observed entries only.
#' @param min_component_support Integer threshold. A missing cell whose
#'   supporting feature has fewer than this many observed entries is labelled
#'   \code{"weakly_estimable"} by \code{\link{get_estimability}} /
#'   \code{\link{get_missing_diagnostics}}.
#' @param rank_candidates Optional non-negative integer vector for
#'   native candidate-rank diagnostics. When \code{NULL}, native automatic rank
#'   selection evaluates \code{0:min(initial_signal_ranks)}.
#' @param n_refits Number of aligned refits for \code{uncertainty =
#'   "bootstrap"} or \code{"mi"}.
#' @param censoring Optional list of censoring metadata, for example
#'   lower-limit vectors by block.
#' @param sensitivity Optional character vector of missingness assumptions to
#'   report in sensitivity diagnostics.
#'
#' @return A list consumed by native missing-data helpers.
#' @export
rajive_missing_control <- function(center = FALSE, scale = FALSE,
                                   normalize = FALSE,
                                   min_component_support = 2L,
                                   rank_candidates = NULL,
                                   n_refits = 5L,
                                   censoring = NULL,
                                   sensitivity = NULL) {
  list(
    center = isTRUE(center),
    scale = isTRUE(scale),
    normalize = isTRUE(normalize),
    min_component_support = as.integer(min_component_support),
    rank_candidates = rank_candidates,
    n_refits = as.integer(n_refits),
    censoring = censoring,
    sensitivity = sensitivity
  )
}

.as_missing_control <- function(control) {
  defaults <- rajive_missing_control()
  if (is.null(control)) return(defaults)
  if (!is.list(control)) {
    cli::cli_abort(
      "`missing_control` must be created by rajive_missing_control().",
      class = "rajiveplus_invalid_input"
    )
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c("Unknown `missing_control` field(s).",
        "x" = "{paste(unknown, collapse = ', ')}."),
      class = "rajiveplus_invalid_input"
    )
  }
  utils::modifyList(defaults, control)
}

.component_decomp_from_matrix <- function(mat, rank = NULL, full = TRUE) {
  if (is.null(rank)) {
    rank <- min(dim(mat))
  }
  rank <- min(as.integer(rank), min(dim(mat)))
  if (rank <= 0L || all(abs(mat) <= .Machine$double.eps)) {
    out <- list(
      u = matrix(0, nrow = nrow(mat), ncol = 0L),
      d = numeric(0L),
      v = matrix(0, nrow = ncol(mat), ncol = 0L),
      full = if (isTRUE(full)) matrix(0, nrow = nrow(mat), ncol = ncol(mat)) else NA,
      rank = 0L
    )
    if (isTRUE(full)) dimnames(out$full) <- dimnames(mat)
    return(out)
  }
  sv <- svd(mat, nu = rank, nv = rank)
  out <- list(
    u = sv$u[, seq_len(rank), drop = FALSE],
    d = sv$d[seq_len(rank)],
    v = sv$v[, seq_len(rank), drop = FALSE],
    full = if (isTRUE(full)) mat else NA,
    rank = rank
  )
  out
}

.component_matrix_from_decomp <- function(decomp, dimnames = NULL) {
  if (is.matrix(decomp$full)) {
    out <- decomp$full
  } else {
    if (is.null(decomp$u) || is.null(decomp$d) || is.null(decomp$v)) {
      cli::cli_abort(
        "Component decomposition does not contain enough SVD information to reconstruct a matrix.",
        class = "rajiveplus_invalid_input"
      )
    }
    out <- svd_reconstruction(decomp)
  }
  if (!is.null(dimnames)) dimnames(out) <- dimnames
  out
}

.fit_component_matrix <- function(fit, k, type, dimnames = NULL) {
  type <- match.arg(type, c("joint", "individual", "noise"))
  if (type == "joint") {
    idx <- 3L * (k - 1L) + 2L
  } else if (type == "individual") {
    idx <- 3L * (k - 1L) + 1L
  } else {
    idx <- 3L * k
  }
  component <- fit$block_decomps[[idx]]
  if (type == "noise" && is.matrix(component)) {
    out <- component
    if (!is.null(dimnames)) dimnames(out) <- dimnames
    return(out)
  }
  .component_matrix_from_decomp(component, dimnames = dimnames)
}

.solve_masked_joint_matrix <- function(x, mask, joint_scores, block_name = NULL) {
  r <- ncol(joint_scores)
  if (r == 0L) {
    out <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    return(out)
  }
  coefs <- matrix(0, nrow = r, ncol = ncol(x))
  ridge <- sqrt(.Machine$double.eps)
  underdetermined <- logical(ncol(x))
  for (j in seq_len(ncol(x))) {
    rows <- mask[, j]
    if (sum(rows) < r) {
      underdetermined[j] <- TRUE
      next
    }
    design <- joint_scores[rows, , drop = FALSE]
    y <- x[rows, j]
    lhs <- crossprod(design) + diag(ridge, r)
    rhs <- crossprod(design, y)
    coefs[, j] <- as.numeric(solve(lhs, rhs))
  }
  if (any(underdetermined)) {
    cli::cli_warn(
      c("Some features have fewer observed entries than the joint rank.",
        "i" = "{if (is.null(block_name)) 'A block' else block_name}: {sum(underdetermined)} feature{?s} had their joint contribution set to zero; that signal is absorbed into the individual component."),
      class = "rajiveplus_underdetermined_joint"
    )
  }
  joint_scores %*% coefs
}

.fit_incomplete_block_decomposition <- function(x, mask, joint_scores,
                                                initial_signal_rank,
                                                full = TRUE,
                                                block_name = NULL) {
  joint_rank <- ncol(joint_scores)
  joint_full <- .solve_masked_joint_matrix(x, mask, joint_scores,
                                           block_name = block_name)
  dimnames(joint_full) <- dimnames(x)
  joint_decomp <- .component_decomp_from_matrix(joint_full,
                                                rank = joint_rank,
                                                full = full)

  residual_seed <- x - joint_full
  residual_seed[!mask] <- NA_real_
  indiv_rank <- max(0L, as.integer(initial_signal_rank) - joint_rank)
  if (indiv_rank > 0L) {
    indiv_svd <- RobRSVD.all(residual_seed, nrank = indiv_rank, weights = mask)
    individual_full <- svd_reconstruction(indiv_svd)
    individual_full[rowSums(mask) == 0L, ] <- 0
    individual_full[, colSums(mask) == 0L] <- 0
    dimnames(individual_full) <- dimnames(x)
    indiv_decomp <- indiv_svd
    indiv_decomp$full <- if (isTRUE(full)) individual_full else NA
    indiv_decomp$rank <- ncol(indiv_svd$u)
  } else {
    individual_full <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    indiv_decomp <- list(
      u = matrix(0, nrow = nrow(x), ncol = 0L),
      d = numeric(0L),
      v = matrix(0, nrow = ncol(x), ncol = 0L),
      full = if (isTRUE(full)) individual_full else NA,
      rank = 0L
    )
  }

  noise <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  noise[mask] <- x[mask] - joint_full[mask] - individual_full[mask]
  list(individual = indiv_decomp, joint = joint_decomp, noise = noise)
}

.attach_native_missing_metadata <- function(fit, normalized, control,
                                            convergence,
                                            preprocess = NULL,
                                            initial_signal_ranks = NULL) {
  fit$missing <- list(
    mask = normalized$mask,
    observed = normalized$observed,
    patterns = normalized$patterns,
    original_blocks = normalized$blocks,
    block_names = normalized$block_names,
    initial_signal_ranks = initial_signal_ranks,
    preprocess = preprocess,
    control = control,
    convergence = convergence
  )
  class(fit) <- unique(c("rajive_incomplete", class(fit)))
  fit$missing$estimability <- .compute_estimability(
    fit, normalized$mask, normalized$patterns
  )
  fit$missing$reconstruction_provenance <- fit$missing$estimability
  fit
}

.Rajive_incomplete <- function(blocks, initial_signal_ranks, joint_rank,
                               mask = NULL,
                               missing_control = rajive_missing_control(),
                               full = TRUE,
                               n_wedin_samples = NA,
                               n_rand_dir_samples = NA,
                               n_perm_samples = NA,
                               num_cores = 1L,
                               seed = NA_integer_,
                               identifiability_norm = "l2",
                               ...) {
  control <- .as_missing_control(missing_control)
  if (!is.na(seed)) set.seed(as.integer(seed))
  .validate_missing_blocks(blocks)
  .validate_native_initial_signal_ranks(initial_signal_ranks, length(blocks))

  auto_rank <- FALSE
  if (is.character(joint_rank)) {
    if (!identical(joint_rank, "native_cv")) {
      cli::cli_abort(
        "`joint_rank` must be numeric or \"native_cv\" in native missing-data mode.",
        class = "rajiveplus_invalid_input"
      )
    }
    auto_rank <- TRUE
  } else if (length(joint_rank) == 1L && is.na(joint_rank)) {
    auto_rank <- TRUE
  }

  if (auto_rank) {
    candidates <- control$rank_candidates
    if (is.null(candidates)) {
      max_candidate <- max(0L, min(as.integer(initial_signal_ranks),
                                   na.rm = TRUE))
      candidates <- seq.int(0L, max_candidate)
    }
    rank_diag <- diagnose_missing_ranks(
      blocks = blocks,
      candidates = candidates,
      initial_signal_ranks = initial_signal_ranks,
      mask = mask,
      missing_control = control,
      full = full,
      num_cores = num_cores,
      seed = NA_integer_,
      identifiability_norm = identifiability_norm
    )
    selected <- rank_diag$joint_rank[which.min(rank_diag$composite_score)]
    fit <- .Rajive_incomplete(
      blocks = blocks,
      initial_signal_ranks = initial_signal_ranks,
      joint_rank = selected,
      mask = mask,
      missing_control = control,
      full = full,
      n_wedin_samples = n_wedin_samples,
      n_rand_dir_samples = n_rand_dir_samples,
      n_perm_samples = n_perm_samples,
      num_cores = num_cores,
      seed = NA_integer_,
      identifiability_norm = identifiability_norm
    )
    fit$missing$rank_diagnostics <- rank_diag
    fit$missing$selected_rank <- selected
    return(fit)
  }
  if (length(joint_rank) != 1L || is.na(joint_rank) ||
      !is.numeric(joint_rank) || joint_rank < 0L) {
    cli::cli_abort(
      "`joint_rank` must be a non-negative integer for native missing-data fits.",
      class = "rajiveplus_invalid_input"
    )
  }
  joint_rank <- as.integer(joint_rank)

  normalized <- .validate_native_missing_inputs(blocks, mask, control)
  mask <- normalized$mask
  any_missing <- any(!unlist(mask, use.names = FALSE))

  if (!any_missing) {
    if (control$center || control$scale || control$normalize) {
      cli::cli_warn(
        c("Native preprocessing was requested but no entries are missing.",
          "i" = "`center`/`scale`/`normalize` apply only when the mask has missing cells; this fit used the unchanged complete-data pipeline."),
        class = "rajiveplus_preprocess_skipped"
      )
    }
    fit <- .Rajive_core(
      blocks = normalized$blocks,
      initial_signal_ranks = initial_signal_ranks,
      full = full,
      n_wedin_samples = n_wedin_samples,
      n_rand_dir_samples = n_rand_dir_samples,
      joint_rank = joint_rank,
      n_perm_samples = n_perm_samples,
      num_cores = num_cores,
      seed = NA_integer_,
      identifiability_norm = identifiability_norm,
      rank_only = FALSE
    )
    block_objective <- vapply(seq_along(normalized$blocks), function(k) {
      joint <- .fit_component_matrix(fit, k, "joint",
                                     dimnames(normalized$blocks[[k]]))
      individual <- .fit_component_matrix(fit, k, "individual",
                                          dimnames(normalized$blocks[[k]]))
      .masked_residual_ss(normalized$blocks[[k]], joint + individual,
                          mask[[k]])
    }, numeric(1L))
    variance <- lapply(seq_along(normalized$blocks), function(k) {
      joint <- .fit_component_matrix(fit, k, "joint",
                                     dimnames(normalized$blocks[[k]]))
      individual <- .fit_component_matrix(fit, k, "individual",
                                          dimnames(normalized$blocks[[k]]))
      .masked_variance_explained(normalized$blocks[[k]], joint, individual,
                                 mask[[k]])
    })
    fit <- .attach_native_missing_metadata(
      fit, normalized, control,
      convergence = list(
        objective = sum(block_objective),
        block_objective = block_objective
      ),
      preprocess = NULL,
      initial_signal_ranks = initial_signal_ranks
    )
    fit$missing$variance_explained <- variance
    return(fit)
  }

  preprocess <- .center_scale_observed(
    normalized$blocks,
    mask,
    center = control$center,
    scale = control$scale,
    normalize = control$normalize
  )
  block_svd <- lapply(seq_along(preprocess$blocks), function(k) {
    get_svd_robustH(preprocess$blocks[[k]],
                    rank = min(initial_signal_ranks[[k]], min(dim(preprocess$blocks[[k]]))),
                    weights = mask[[k]])
  })

  signal_scores <- lapply(seq_along(block_svd), function(k) {
    rank_k <- min(initial_signal_ranks[[k]], ncol(block_svd[[k]]$u))
    block_svd[[k]]$u[, seq_len(rank_k), drop = FALSE]
  })
  signal_matrix <- do.call(cbind, signal_scores)
  rank_for_signal <- min(joint_rank, min(dim(signal_matrix)))
  if (rank_for_signal == 0L) {
    joint_scores <- matrix(0, nrow = nrow(signal_matrix), ncol = 0L)
    obs_svals <- numeric(0L)
  } else {
    signal_svd <- get_svd_robustH(signal_matrix, rank = rank_for_signal)
    joint_scores <- signal_svd$u[, seq_len(rank_for_signal), drop = FALSE]
    obs_svals <- signal_svd$d
  }
  if (ncol(joint_scores) < joint_rank) {
    cli::cli_warn(
      c("Requested `joint_rank` exceeds the recoverable joint dimension.",
        "i" = "Reduced from {joint_rank} to {ncol(joint_scores)} (rank of the concatenated signal scores)."),
      class = "rajiveplus_joint_rank_truncated"
    )
    joint_rank <- ncol(joint_scores)
  }

  block_decomps <- vector("list", length(preprocess$blocks) * 3L)
  objectives <- numeric(length(preprocess$blocks))
  variance <- vector("list", length(preprocess$blocks))
  for (k in seq_along(preprocess$blocks)) {
    decomp <- .fit_incomplete_block_decomposition(
      preprocess$blocks[[k]], mask[[k]], joint_scores,
      initial_signal_rank = initial_signal_ranks[[k]],
      full = full,
      block_name = preprocess$block_names[[k]]
    )
    block_decomps[[3L * (k - 1L) + 1L]] <- decomp$individual
    block_decomps[[3L * (k - 1L) + 2L]] <- decomp$joint
    block_decomps[[3L * k]] <- decomp$noise
    joint_full <- .component_matrix_from_decomp(
      decomp$joint, dimnames(preprocess$blocks[[k]])
    )
    individual_full <- .component_matrix_from_decomp(
      decomp$individual, dimnames(preprocess$blocks[[k]])
    )
    fitted <- joint_full + individual_full
    objectives[[k]] <- .masked_residual_ss(preprocess$blocks[[k]], fitted, mask[[k]])
    variance[[k]] <- .masked_variance_explained(preprocess$blocks[[k]],
                                                joint_full,
                                                individual_full,
                                                mask[[k]])
  }

  fit <- list(
    block_decomps = block_decomps,
    joint_scores = joint_scores,
    joint_rank = joint_rank,
    joint_rank_sel = list(
      obs_svals = obs_svals,
      joint_rank_estimate = joint_rank,
      overall_sv_sq_threshold = NA_real_,
      identifiability_norm = identifiability_norm,
      native_missing = TRUE
    )
  )
  class(fit) <- "rajive"

  convergence <- list(
    objective = sum(objectives),
    block_objective = objectives
  )
  fit <- .attach_native_missing_metadata(fit, normalized, control,
                                         convergence, preprocess,
                                         initial_signal_ranks = initial_signal_ranks)
  fit$missing$variance_explained <- variance
  fit
}

.compute_estimability <- function(fit, mask, patterns, diagnostics = NULL) {
  block_names <- fit$missing$block_names
  if (is.null(block_names)) block_names <- names(mask)
  if (is.null(block_names) || any(!nzchar(block_names))) {
    block_names <- paste0("block", seq_along(mask))
  }
  min_support <- fit$missing$control$min_component_support
  if (is.null(min_support) || length(min_support) == 0L ||
      !is.finite(min_support[[1L]])) {
    min_support <- 1L
  } else {
    min_support <- as.integer(min_support[[1L]])
  }
  sample_support <- patterns$sample_summary$observed_cells > 0L
  rows <- lapply(seq_along(mask), function(k) {
    obs <- mask[[k]]
    whole_row <- rowSums(obs) == 0L
    partial_row <- rowSums(obs) > 0L & rowSums(obs) < ncol(obs)
    feature_obs <- colSums(obs)
    feature_support <- feature_obs > 0L
    feature_strong <- feature_obs >= min_support
    grid <- expand.grid(
      row = seq_len(nrow(obs)),
      col = seq_len(ncol(obs)),
      component = c("joint", "individual", "joint_individual"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    observed <- obs[cbind(grid$row, grid$col)]
    label <- rep("observed", nrow(grid))
    missing_cell <- !observed

    joint_only <- grid$component == "joint"
    indiv_only <- grid$component == "individual"
    joint_indiv <- grid$component == "joint_individual"
    whole_row_cell <- whole_row[grid$row]
    partial_row_cell <- partial_row[grid$row]
    feature_supported <- feature_support[grid$col]
    feature_strong_cell <- feature_strong[grid$col]
    cross_supported <- sample_support[grid$row] & fit$joint_rank > 0L

    label[missing_cell & whole_row_cell & joint_only & cross_supported] <-
      "cross_block_joint_estimable"
    label[missing_cell & whole_row_cell & joint_indiv & cross_supported] <-
      "cross_block_joint_estimable"
    label[missing_cell & whole_row_cell & indiv_only] <- "not_identifiable"
    label[missing_cell & whole_row_cell & !cross_supported & !indiv_only] <-
      "not_identifiable"

    within <- missing_cell & !whole_row_cell & partial_row_cell &
      feature_strong_cell
    label[within] <- "within_block_estimable"

    weak <- missing_cell & !whole_row_cell & partial_row_cell &
      feature_supported & !feature_strong_cell
    label[weak] <- "weakly_estimable"
    label[missing_cell & !whole_row_cell & !partial_row_cell] <- "not_identifiable"

    data.frame(
      block = block_names[[k]],
      row = grid$row,
      col = grid$col,
      sample_id = {
        rn <- rownames(obs)
        if (is.null(rn)) as.character(grid$row) else rn[grid$row]
      },
      feature_id = {
        cn <- colnames(obs)
        if (is.null(cn)) as.character(grid$col) else cn[grid$col]
      },
      component = grid$component,
      observed = observed,
      label = label,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.resolve_missing_block_filter <- function(block, block_names) {
  if (is.null(block)) return(seq_along(block_names))
  if (is.character(block)) {
    idx <- match(block, block_names)
  } else {
    idx <- as.integer(block)
  }
  if (anyNA(idx) || any(idx < 1L | idx > length(block_names))) {
    cli::cli_abort(
      "`block` must identify an existing block.",
      class = "rajiveplus_invalid_input"
    )
  }
  idx
}

#' Get native missing-data estimability labels
#'
#' @param fit A \code{rajive_incomplete} object returned by
#'   \code{\link{Rajive}} with \code{missing = "native"}.
#' @param block Optional block index or name.
#' @param component Optional component label: \code{"joint"},
#'   \code{"individual"}, or \code{"joint_individual"}.
#'
#' @return A data frame with one row per block/sample/feature/component label.
#' @export
get_estimability <- function(fit, block = NULL, component = NULL) {
  if (is.null(fit$missing)) {
    cli::cli_abort(
      "`fit` does not contain native missing-data metadata.",
      class = "rajiveplus_invalid_input"
    )
  }
  est <- fit$missing$estimability
  if (is.null(est)) {
    est <- .compute_estimability(fit, fit$missing$mask, fit$missing$patterns)
  }
  block_names <- fit$missing$block_names
  idx <- .resolve_missing_block_filter(block, block_names)
  est <- est[est$block %in% block_names[idx], , drop = FALSE]
  if (!is.null(component)) {
    est <- est[est$component %in% component, , drop = FALSE]
  }
  rownames(est) <- NULL
  est
}

.standardized_original_blocks <- function(fit) {
  if (!is.null(fit$missing$preprocess)) {
    return(fit$missing$preprocess$blocks)
  }
  fit$missing$original_blocks
}

#' Reconstruct blocks from a native missing-data fit
#'
#' @param fit A \code{rajive_incomplete} object.
#' @param type \code{"joint"} or \code{"joint_individual"}.
#' @param scale \code{"original"} or \code{"standardized"}.
#' @param observed \code{"fitted"} returns fitted values at observed entries;
#'   \code{"observed"} restores measured values at observed entries.
#'
#' @return A named list of reconstructed matrices. The returned list has a
#'   \code{"reconstruction_provenance"} attribute containing matching
#'   estimability labels.
#' @export
get_reconstructed_blocks <- function(fit,
                                     type = c("joint", "joint_individual"),
                                     scale = c("original", "standardized"),
                                     observed = c("fitted", "observed")) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  observed <- match.arg(observed)
  if (is.null(fit$missing)) {
    cli::cli_abort(
      "`fit` does not contain native missing-data metadata.",
      class = "rajiveplus_invalid_input"
    )
  }
  K <- length(fit$missing$mask)
  block_names <- fit$missing$block_names
  out <- vector("list", K)
  names(out) <- block_names
  for (k in seq_len(K)) {
    joint <- .fit_component_matrix(fit, k, "joint",
                                   dimnames(fit$missing$mask[[k]]))
    if (type == "joint") {
      mat <- joint
    } else {
      individual <- .fit_component_matrix(fit, k, "individual",
                                          dimnames(fit$missing$mask[[k]]))
      mat <- joint + individual
      whole_missing_rows <- rowSums(fit$missing$mask[[k]]) == 0L
      if (any(whole_missing_rows)) {
        mat[whole_missing_rows, ] <- joint[whole_missing_rows, ]
      }
    }
    dimnames(mat) <- dimnames(fit$missing$mask[[k]])
    out[[k]] <- mat
  }

  if (scale == "original" && !is.null(fit$missing$preprocess)) {
    out <- .backtransform_reconstruction(out, fit$missing$preprocess)
  }

  if (observed == "observed") {
    observed_blocks <- if (scale == "standardized") {
      .standardized_original_blocks(fit)
    } else {
      fit$missing$original_blocks
    }
    for (k in seq_len(K)) {
      obs <- fit$missing$mask[[k]]
      out[[k]][obs] <- observed_blocks[[k]][obs]
    }
  }

  attr(out, "reconstruction_provenance") <- get_estimability(
    fit, component = type
  )
  out
}

#' Get native missingness summaries
#'
#' @param fit A \code{rajive_incomplete} object.
#' @return A list of block, sample, feature, cell, and count summaries.
#' @export
get_missingness_info <- function(fit) {
  if (is.null(fit$missing)) {
    cli::cli_abort(
      "`fit` does not contain native missing-data metadata.",
      class = "rajiveplus_invalid_input"
    )
  }
  fit$missing$patterns
}

.missing_variance_table <- function(fit) {
  variance <- fit$missing$variance_explained
  block_names <- fit$missing$block_names
  if (is.null(variance)) {
    variance <- lapply(seq_along(fit$missing$mask), function(k) {
      .masked_variance_explained(
        .standardized_original_blocks(fit)[[k]],
        .fit_component_matrix(fit, k, "joint",
                              dimnames(fit$missing$mask[[k]])),
        .fit_component_matrix(fit, k, "individual",
                              dimnames(fit$missing$mask[[k]])),
        fit$missing$mask[[k]]
      )
    })
  }
  out <- do.call(rbind, lapply(seq_along(variance), function(k) {
    data.frame(
      block = block_names[[k]],
      as.list(variance[[k]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

.joint_support_table <- function(fit) {
  mask <- fit$missing$mask
  block_names <- fit$missing$block_names
  n <- nrow(mask[[1L]])
  r <- max(1L, fit$joint_rank)
  rows <- lapply(seq_len(n), function(i) {
    observed_by_block <- vapply(mask, function(m) sum(m[i, ]), integer(1L))
    data.frame(
      sample = i,
      component = seq_len(r),
      observed_blocks = sum(observed_by_block > 0L),
      observed_cells = sum(observed_by_block),
      weak_support = sum(observed_by_block > 0L) < 2L,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$block_count <- length(block_names)
  rownames(out) <- NULL
  out
}

#' Get native missing-data diagnostics
#'
#' @param fit A \code{rajive_incomplete} object.
#' @param warn_weak Logical. If \code{TRUE}, emit a classed warning when any
#'   quantity is labelled \code{"weakly_estimable"}.
#'
#' @return A list containing missingness, estimability, variance, support,
#'   orthogonality, convergence, assumption, censoring, and sensitivity
#'   diagnostics.
#' @export
get_missing_diagnostics <- function(fit, warn_weak = FALSE) {
  if (is.null(fit$missing)) {
    cli::cli_abort(
      "`fit` does not contain native missing-data metadata.",
      class = "rajiveplus_invalid_input"
    )
  }
  est <- get_estimability(fit)
  if (isTRUE(warn_weak) && any(est$label == "weakly_estimable")) {
    cli::cli_warn(
      "Some native missing-data quantities are weakly estimable.",
      class = "rajiveplus_weak_estimability"
    )
  }
  list(
    missingness = get_missingness_info(fit),
    estimability = est,
    variance_explained = .missing_variance_table(fit),
    joint_support = .joint_support_table(fit),
    orthogonality = do.call(rbind, lapply(seq_along(fit$missing$mask), function(k) {
      dn <- dimnames(fit$missing$mask[[k]])
      joint <- .fit_component_matrix(fit, k, "joint", dn)
      individual <- .fit_component_matrix(fit, k, "individual", dn)
      o <- .masked_orthogonality(joint, individual, fit$missing$mask[[k]])
      data.frame(block = fit$missing$block_names[[k]],
                 max_abs_crossprod = o$max_abs_crossprod,
                 observed_fraction = o$observed_fraction,
                 stringsAsFactors = FALSE)
    })),
    convergence = fit$missing$convergence,
    assumptions = paste(
      "Native observed-entry fitting treats missingness as ignorable for the",
      "scientific estimand, broadly MCAR/MAR-style; generic native mode is",
      "not a general MNAR solution."
    ),
    censoring = .missing_censoring_summary(fit),
    sensitivity = .missing_sensitivity_summary(fit)
  )
}

#' Plot native missingness patterns
#'
#' @param fit A \code{rajive_incomplete} object.
#' @return A \code{ggplot} heatmap of observed versus missing cells by block.
#' @export
plot_missingness <- function(fit) {
  info <- get_missingness_info(fit)
  cells <- info$cells
  ggplot2::ggplot(cells, ggplot2::aes(x = .data$col, y = .data$row,
                                      fill = .data$observed)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~ block, scales = "free_x") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#2c7bb6",
                                          "FALSE" = "#d7191c")) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "Feature", y = "Sample", fill = "Observed") +
    ggplot2::theme_bw()
}

#' @export
print.rajive_incomplete <- function(x, ...) {
  NextMethod("print")
  info <- get_missingness_info(x)
  cat(sprintf("  Missing cells    : %d\n", info$counts[["missing_cells"]]))
  cat(sprintf("  Observed fraction: %.3f\n",
              sum(info$block_summary$observed_cells) /
                sum(info$block_summary$observed_cells + info$block_summary$missing_cells)))
  invisible(x)
}

#' @export
summary.rajive_incomplete <- function(object, ...) {
  ranks <- NextMethod("summary")
  diag <- get_missing_diagnostics(object)
  cat("\nNative missing-data diagnostics\n")
  print(diag$variance_explained, row.names = FALSE)
  invisible(list(ranks = ranks, missing = diag))
}

#' Diagnose candidate joint ranks for native missing-data fits
#'
#' @param blocks List of numeric data matrices.
#' @param candidates Non-negative candidate joint ranks.
#' @param initial_signal_ranks Initial signal ranks, one per block.
#' @param mask Optional observation mask.
#' @param missing_control Control list from
#'   \code{\link{rajive_missing_control}}.
#' @param full,num_cores,seed,identifiability_norm Arguments forwarded to the
#'   native fixed-rank fit.
#' @param ... Reserved for future diagnostics.
#'
#' @details Native automatic rank selection picks the candidate minimising
#'   \code{composite_score} = observed-entry \code{prediction_error} plus a
#'   parsimony penalty of 0.001 per joint component.
#'
#' @return A data frame with columns \code{joint_rank},
#'   \code{prediction_error}, \code{weak_support_rate},
#'   \code{missing_fraction}, and \code{composite_score}.
#' @export
diagnose_missing_ranks <- function(blocks,
                                   candidates,
                                   initial_signal_ranks,
                                   mask = NULL,
                                   missing_control = rajive_missing_control(),
                                   full = TRUE,
                                   num_cores = 1L,
                                   seed = NA_integer_,
                                   identifiability_norm = "l2",
                                   ...) {
  control <- .as_missing_control(missing_control)
  normalized <- .validate_native_missing_inputs(blocks, mask, control)
  .validate_native_initial_signal_ranks(initial_signal_ranks,
                                        length(normalized$blocks))
  candidates <- sort(unique(as.integer(candidates)))
  if (length(candidates) == 0L || anyNA(candidates) || any(candidates < 0L)) {
    cli::cli_abort(
      "`candidates` must contain non-negative integer ranks.",
      class = "rajiveplus_invalid_input"
    )
  }
  max_rank <- min(initial_signal_ranks)
  candidates <- candidates[candidates <= max_rank]
  if (length(candidates) == 0L) {
    cli::cli_abort(
      "No candidate rank is compatible with `initial_signal_ranks`.",
      class = "rajiveplus_invalid_input"
    )
  }
  if (!is.na(seed)) set.seed(as.integer(seed))

  missing_fraction <- mean(!unlist(normalized$mask, use.names = FALSE))
  diagnostic_blocks <- normalized$blocks
  if (missing_fraction > 0) {
    diagnostic_blocks <- .center_scale_observed(
      normalized$blocks,
      normalized$mask,
      center = control$center,
      scale = control$scale,
      normalize = control$normalize
    )$blocks
  }
  observed_ss <- sum(vapply(seq_along(diagnostic_blocks), function(k) {
    sum(diagnostic_blocks[[k]][normalized$mask[[k]]]^2)
  }, numeric(1L)))
  if (!is.finite(observed_ss) || observed_ss <= .Machine$double.eps) observed_ss <- 1

  rows <- lapply(candidates, function(rank) {
    fit <- .Rajive_incomplete(
      blocks = normalized$blocks,
      initial_signal_ranks = initial_signal_ranks,
      joint_rank = rank,
      mask = normalized$mask,
      missing_control = control,
      full = full,
      num_cores = num_cores,
      seed = NA_integer_,
      identifiability_norm = identifiability_norm
    )
    pred <- fit$missing$convergence$objective / observed_ss
    support <- .joint_support_table(fit)
    weak_rate <- mean(support$weak_support)
    data.frame(
      joint_rank = rank,
      prediction_error = pred,
      weak_support_rate = weak_rate,
      missing_fraction = missing_fraction,
      composite_score = pred + 0.001 * rank,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.compute_missing_uncertainty <- function(fit, method = c("bootstrap", "mi")) {
  method <- match.arg(method)
  control <- fit$missing$control
  n_refits <- max(1L, as.integer(control$n_refits))
  base_blocks <- fit$missing$original_blocks
  mask <- fit$missing$mask
  initial_signal_ranks <- fit$missing$initial_signal_ranks
  joint_rank <- fit$joint_rank
  refits <- vector("list", n_refits)
  recon <- get_reconstructed_blocks(fit, type = "joint_individual",
                                    scale = "original")

  for (b in seq_len(n_refits)) {
    refit_control <- control
    if (method == "bootstrap") {
      perturbed <- lapply(seq_along(base_blocks), function(k) {
        x <- base_blocks[[k]]
        obs <- mask[[k]]
        sigma <- stats::sd(x[obs] - recon[[k]][obs], na.rm = TRUE)
        if (!is.finite(sigma) || sigma <= 0) sigma <- 1e-6
        x[obs] <- x[obs] + stats::rnorm(sum(obs), sd = sigma)
        x
      })
      names(perturbed) <- names(base_blocks)
      refit_mask <- mask
    } else {
      perturbed <- lapply(seq_along(base_blocks), function(k) {
        x <- base_blocks[[k]]
        obs <- mask[[k]]
        sigma <- stats::sd(x[obs] - recon[[k]][obs], na.rm = TRUE)
        if (!is.finite(sigma) || sigma <= 0) sigma <- 1e-6
        x[!obs] <- recon[[k]][!obs] + stats::rnorm(sum(!obs), sd = sigma)
        x
      })
      names(perturbed) <- names(base_blocks)
      refit_mask <- lapply(perturbed, function(x) matrix(TRUE, nrow(x), ncol(x)))
      refit_control <- utils::modifyList(
        control,
        list(center = FALSE, scale = FALSE, normalize = FALSE)
      )
    }
    refit <- .Rajive_incomplete(
      perturbed,
      initial_signal_ranks = initial_signal_ranks,
      joint_rank = joint_rank,
      mask = refit_mask,
      missing_control = refit_control,
      full = TRUE
    )
    refits[[b]] <- .align_missing_refit(fit, refit)
  }

  score_summary <- .summarize_missing_refit_scores(fit, refits)
  list(
    method = method,
    refits = refits,
    score_summary = score_summary,
    adapter_note = "Use get_missing_uncertainty() summaries as a native-missing adapter for rajive_ci-style reporting; raw unaligned refits are not pooled."
  )
}

.align_missing_refit <- function(reference, refit) {
  r <- min(ncol(reference$joint_scores), ncol(refit$joint_scores))
  if (r > 0L) {
    for (j in seq_len(r)) {
      cc <- suppressWarnings(stats::cor(reference$joint_scores[, j],
                                        refit$joint_scores[, j]))
      if (is.finite(cc) && cc < 0) {
        refit$joint_scores[, j] <- -refit$joint_scores[, j]
        K <- length(refit$missing$mask)
        for (k in seq_len(K)) {
          idx <- 3L * (k - 1L) + 2L
          if (ncol(refit$block_decomps[[idx]]$u) >= j) {
            refit$block_decomps[[idx]]$u[, j] <- -refit$block_decomps[[idx]]$u[, j]
          }
        }
      }
    }
  }
  refit
}

.summarize_missing_refit_scores <- function(fit, refits) {
  r <- ncol(fit$joint_scores)
  if (r == 0L) {
    return(data.frame(
      sample = integer(0),
      component = integer(0),
      sd = numeric(0),
      missing_blocks = integer(0)
    ))
  }
  rows <- list()
  mask <- fit$missing$mask
  missing_blocks <- Reduce(`+`, lapply(mask, function(m) rowSums(m) == 0L))
  for (j in seq_len(r)) {
    score_mat <- do.call(cbind, lapply(refits, function(x) x$joint_scores[, j]))
    rows[[j]] <- data.frame(
      sample = seq_len(nrow(score_mat)),
      component = j,
      sd = apply(score_mat, 1L, stats::sd),
      missing_blocks = as.integer(missing_blocks),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Get native missing-data uncertainty summaries
#'
#' @param fit A \code{rajive_incomplete} object fit with
#'   \code{uncertainty = "bootstrap"} or \code{"mi"}.
#' @return A list containing aligned refits and score summaries.
#' @export
get_missing_uncertainty <- function(fit) {
  if (is.null(fit$missing) || is.null(fit$missing$uncertainty)) {
    cli::cli_abort(
      "`fit` does not contain native missing-data uncertainty results.",
      class = "rajiveplus_invalid_input"
    )
  }
  fit$missing$uncertainty
}

.missing_censoring_summary <- function(fit) {
  censoring <- fit$missing$control$censoring
  if (is.null(censoring)) {
    return(list(enabled = FALSE, blocks = character(0), note = NULL))
  }
  lower <- censoring$lower_limit
  blocks <- if (is.list(lower)) names(lower) else fit$missing$block_names
  blocks <- blocks[nzchar(blocks)]
  list(
    enabled = TRUE,
    blocks = blocks,
    note = "Censoring metadata is used for diagnostics/sensitivity reporting; generic native fitting remains observed-entry based."
  )
}

.missing_sensitivity_summary <- function(fit) {
  assumptions <- fit$missing$control$sensitivity
  if (is.null(assumptions)) {
    return(list(assumptions = character(0), implemented = FALSE, note = NULL))
  }
  list(
    assumptions = as.character(assumptions),
    implemented = FALSE,
    note = paste(
      "Sensitivity analysis under alternative missingness assumptions is not",
      "yet implemented; the requested assumptions are recorded for reference",
      "only and do not alter the fit."
    )
  )
}
