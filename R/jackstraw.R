# Jackstraw significance testing for RaJIVE joint loadings.
#
# Implements the permutation-based jackstraw test described in:
#   Yang X, Hoadley KA, Hannig J, Marron JS (2021). Statistical inference for
#   data integration. arXiv:2109.12272.
#
# The test identifies which variables (features) in each data block have
# statistically significantly nonzero joint loadings after a RaJIVE
# decomposition, adapting the jackstraw idea of Chung & Storey (2015) to
# the JIVE setting.

# ---------------------------------------------------------------------------
# Part 1: Core statistical helpers (internal, not exported)
# ---------------------------------------------------------------------------

#' Vectorized OLS F-statistics for simple linear regression
#'
#' For every row of \code{Y_t} (each feature), compute the F-statistic for
#' \code{feature ~ x + 1} in a single matrix operation.  Constant rows
#' (zero variance) return \code{NA}.
#'
#' @param Y_t \code{d x n} numeric matrix; features are rows, observations
#'   are columns.
#' @param x length-\code{n} numeric vector of predictor values (joint scores
#'   for one component).
#'
#' @return length-\code{d} numeric vector of F-statistics (\code{NA} for
#'   constant rows).
#'
#' @keywords internal
ols_f_stat_matrix <- function(Y_t, x) {
  n   <- length(x)
  x_c <- x - mean(x)
  # center each feature row
  row_means <- rowMeans(Y_t)                      # d
  Y_c       <- Y_t - row_means                    # d x n

  xTx <- sum(x_c^2)

  # beta1 for all features at once: (Y_c %*% x_c) / xTx
  beta1 <- as.numeric(Y_c %*% x_c) / xTx         # d

  # residuals: Y_c - beta1 * x_c'  (broadcasting via outer)
  E   <- Y_c - tcrossprod(beta1, x_c)             # d x n

  sse1 <- rowSums(E^2)                            # d
  ss0  <- rowSums(Y_c^2)                          # d

  # Features with zero variance: sse1 == ss0, mark as NA
  F_stat <- (ss0 - sse1) / (sse1 / (n - 2))

  # Constant features: use a tolerance rather than exact equality to catch
  # near-zero variances caused by floating-point arithmetic.
  is_constant <- (rowVars_fast(Y_t) < .Machine$double.eps)
  F_stat[is_constant] <- NA_real_

  F_stat
}

#' Fast row variances (internal helper)
#'
#' @param M numeric matrix
#' @return numeric vector of row variances
#' @keywords internal
rowVars_fast <- function(M) {
  mu  <- rowMeans(M)
  rowSums((M - mu)^2) / (ncol(M) - 1)
}


#' Compute empirical p-values from observed and null F-statistics
#'
#' For feature \code{i}, the p-value is the proportion of its \code{n_null}
#' null F-statistics that exceed the observed value.  \code{NA} observed
#' statistics (constant features) receive p-value = 1.
#'
#' @param f_obs length-\code{d} numeric vector of observed F-statistics.
#' @param f_null \code{d x n_null} numeric matrix of null F-statistics.
#'
#' @return length-\code{d} numeric vector of empirical p-values in [0, 1].
#'
#' @keywords internal
compute_empirical_pvalues <- function(f_obs, f_null) {
  d      <- length(f_obs)
  p_vals <- numeric(d)

  for (i in seq_len(d)) {
    if (is.na(f_obs[i])) {
      p_vals[i] <- 1
    } else {
      nulls      <- f_null[i, ]
      nulls      <- nulls[!is.na(nulls)]
      if (length(nulls) == 0L) {
        p_vals[i] <- 1
      } else {
        p_vals[i] <- mean(f_obs[i] < nulls)
      }
    }
  }
  p_vals
}


#' Generate null F-statistics via jackstraw permutation
#'
#' For each of the \code{d} features, sample \code{n_null} feature rows from
#' \code{X_t} (with replacement), permute the sampled rows' values (breaking
#' any association with \code{joint_comp_scores}), and compute their
#' F-statistics.
#'
#' @param X_t \code{d x n} numeric matrix (data block, transposed).
#' @param joint_comp_scores length-\code{n} numeric vector of joint scores for
#'   one component.
#' @param n_null positive integer; number of null F-statistics per feature.
#'
#' @return \code{d x n_null} numeric matrix of null F-statistics.
#'
#' @keywords internal
generate_null_f_stats <- function(X_t, joint_comp_scores, n_null) {
  d <- nrow(X_t)
  n <- ncol(X_t)

  # Pre-sample feature indices for all (d * n_null) draws at once
  sampled_idx <- matrix(
    sample.int(d, size = d * n_null, replace = TRUE),
    nrow = d, ncol = n_null
  )

  null_f <- matrix(NA_real_, nrow = d, ncol = n_null)

  for (s in seq_len(n_null)) {
    # Extract sampled rows (d x n).  Instead of independently permuting each
    # sampled feature row (which requires d separate permutations), we permute
    # the score vector once.  Permuting scores breaks the feature-score
    # association in the same way and yields the same null distribution:
    #   F(feature_perm, scores) =_d F(feature, scores_perm)  under H0.
    rows        <- X_t[sampled_idx[, s], , drop = FALSE]   # d x n
    scores_perm <- sample(joint_comp_scores)                # single permutation
    null_f[, s] <- ols_f_stat_matrix(rows, scores_perm)
  }

  null_f
}


# ---------------------------------------------------------------------------
# Part 2: Main exported function
# ---------------------------------------------------------------------------

#' Jackstraw significance testing for RaJIVE joint loadings
#'
#' Applies a permutation-based jackstraw test to determine which features in
#' each data block have statistically significantly nonzero joint loadings
#' from a \code{\link{Rajive}} decomposition.
#'
#' For each data block \eqn{k} and each joint component \eqn{j}, the observed
#' F-statistic for the regression \emph{feature ~ joint_score_j + 1} is
#' compared to a null distribution generated by permuting randomly sampled
#' feature values, thereby breaking the association with the joint scores.
#' Empirical p-values are computed and optionally corrected for multiple
#' testing.
#'
#' @param ajive_output List returned by \code{\link{Rajive}}.
#' @param blocks List of data matrices (same list passed to
#'   \code{\link{Rajive}}).
#' @param alpha Numeric scalar; desired significance level.  Default
#'   \code{0.05}.
#' @param n_null Positive integer; number of null F-statistics generated per
#'   feature per joint component.  Larger values give more stable p-values
#'   at the cost of computation time.  Default \code{10}; recommended
#'   \code{50}--\code{100} for publication-quality results.
#' @param correction Character string controlling multiple-testing correction.
#'   One of \code{"bonferroni"} (default), \code{"BH"} (Benjamini--Hochberg
#'   FDR), or \code{"none"}.  \code{"bonferroni"} divides \code{alpha} by
#'   \eqn{d_k \times \text{joint\_rank}} for each block, matching the
#'   original Python implementation.
#'
#' @return An object of class \code{"jackstraw_rajive"}: a named list with one
#'   element per block (\code{block1}, \code{block2}, \ldots).  Each element
#'   is itself a list with one element per joint component
#'   (\code{comp1}, \code{comp2}, \ldots) containing:
#'   \describe{
#'     \item{\code{f_obs}}{length-\eqn{d_k} numeric vector of observed
#'       F-statistics.}
#'     \item{\code{f_null}}{\eqn{d_k \times} \code{n_null} matrix of null
#'       F-statistics.}
#'     \item{\code{p_values}}{Empirical p-values (length \eqn{d_k}).}
#'     \item{\code{p_adj}}{Multiple-testing adjusted p-values (length
#'       \eqn{d_k}).}
#'     \item{\code{significant}}{Named logical vector (length \eqn{d_k})
#'       indicating significance.}
#'     \item{\code{significant_vars}}{Integer indices (or column names when
#'       available) of significant features.}
#'   }
#'   The object also carries attributes \code{alpha}, \code{correction},
#'   \code{joint_rank}, and \code{n_blocks}.
#'
#' @references
#' Yang X, Hoadley KA, Hannig J, Marron JS (2021). Statistical inference for
#' data integration. \emph{arXiv:2109.12272}.
#'
#' Chung NC, Storey JD (2015). Statistical significance of variables driving
#' systematic variation in high-dimensional data.
#' \emph{Bioinformatics}, 31(4):545--554.
#'
#' @seealso \code{\link{Rajive}}, \code{\link{plot_jackstraw}},
#'   \code{\link{get_significant_vars}}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n   <- 50
#' pks <- c(100, 80)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = n,
#'                       pks = pks, dist.type = 1)
#' data.ajive           <- Y$sim_data
#' initial_signal_ranks <- c(5, 4)
#' ajive_result <- Rajive(data.ajive, initial_signal_ranks)
#' js <- jackstraw_rajive(ajive_result, data.ajive, alpha = 0.05, n_null = 10)
#' print(js)
#' summary(js)
#' }
#'
#' @export
jackstraw_rajive <- function(ajive_output, blocks,
                             alpha      = 0.05,
                             n_null     = 10,
                             correction = c("bonferroni", "BH", "none")) {

  correction <- match.arg(correction)

  if (!is.list(blocks) || length(blocks) == 0L) {
    stop("'blocks' must be a non-empty list of data matrices.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(n_null) || length(n_null) != 1L || n_null < 1L) {
    stop("'n_null' must be a positive integer.", call. = FALSE)
  }
  n_null <- as.integer(n_null)

  joint_scores <- ajive_output[["joint_scores"]]
  if (is.null(joint_scores)) {
    stop("'ajive_output' does not contain 'joint_scores'. ",
         "Make sure it is the output of Rajive().", call. = FALSE)
  }
  joint_rank <- ncol(joint_scores)
  K          <- length(blocks)

  if (joint_rank == 0L) {
    stop("The estimated joint rank is 0; jackstraw testing is not applicable.",
         call. = FALSE)
  }

  results <- vector("list", K)
  names(results) <- paste0("block", seq_len(K))

  for (k in seq_len(K)) {
    X   <- blocks[[k]]
    X_t <- t(X)          # d x n
    d_k <- ncol(X)

    # Bonferroni: divide alpha by total number of tests in this block
    alpha_bonf <- if (correction == "bonferroni") {
      alpha / (d_k * joint_rank)
    } else {
      alpha
    }

    block_results <- vector("list", joint_rank)
    names(block_results) <- paste0("comp", seq_len(joint_rank))

    for (j in seq_len(joint_rank)) {
      scores_j <- joint_scores[, j]

      f_obs  <- ols_f_stat_matrix(X_t, scores_j)
      f_null <- generate_null_f_stats(X_t, scores_j, n_null)
      p_vals <- compute_empirical_pvalues(f_obs, f_null)

      # Multiple-testing correction
      if (correction == "BH") {
        p_adj <- stats::p.adjust(p_vals, method = "BH")
        sig   <- p_adj <= alpha
      } else if (correction == "bonferroni") {
        p_adj <- p_vals
        sig   <- p_vals <= alpha_bonf
      } else {
        # "none"
        p_adj <- p_vals
        sig   <- p_vals <= alpha
      }

      # Attach variable names where available
      feat_names <- colnames(X)
      if (!is.null(feat_names)) {
        names(f_obs)  <- feat_names
        names(p_vals) <- feat_names
        names(p_adj)  <- feat_names
        names(sig)    <- feat_names
        sig_vars      <- feat_names[sig & !is.na(sig)]
      } else {
        sig_vars <- which(sig & !is.na(sig))
      }

      block_results[[j]] <- list(
        f_obs           = f_obs,
        f_null          = f_null,
        p_values        = p_vals,
        p_adj           = p_adj,
        significant     = sig,
        significant_vars = sig_vars
      )
    }

    results[[k]] <- block_results
  }

  attr(results, "alpha")      <- alpha
  attr(results, "correction") <- correction
  attr(results, "joint_rank") <- joint_rank
  attr(results, "n_blocks")   <- K

  class(results) <- "jackstraw_rajive"
  results
}


# ---------------------------------------------------------------------------
# Part 3: S3 methods
# ---------------------------------------------------------------------------

#' Print method for jackstraw_rajive objects
#'
#' Displays a concise summary table of significance results for each block and
#' joint component.
#'
#' @param x An object of class \code{"jackstraw_rajive"}.
#' @param ... Ignored.
#'
#' @return \code{x} invisibly.
#'
#' @export
print.jackstraw_rajive <- function(x, ...) {
  alpha      <- attr(x, "alpha")
  correction <- attr(x, "correction")
  joint_rank <- attr(x, "joint_rank")
  K          <- attr(x, "n_blocks")

  cat("JIVE Jackstraw Significance Test\n")
  cat(sprintf("  Joint rank: %d   Alpha: %g   Correction: %s\n\n",
              joint_rank, alpha, correction))
  cat(sprintf("  %-10s %-12s %-14s %-14s\n",
              "Block", "Component", "N features", "N significant"))
  cat("  ", strrep("-", 52), "\n", sep = "")

  for (k in seq_len(K)) {
    block_nm <- names(x)[k]
    for (j in seq_len(joint_rank)) {
      comp   <- x[[k]][[j]]
      n_feat <- length(comp$f_obs)
      n_sig  <- sum(comp$significant, na.rm = TRUE)
      cat(sprintf("  %-10s %-12s %-14d %-14d\n",
                  block_nm, paste0("comp", j), n_feat, n_sig))
    }
  }
  invisible(x)
}


#' Summary method for jackstraw_rajive objects
#'
#' Returns (and prints) a \code{data.frame} with one row per block/component
#' combination.
#'
#' @param object An object of class \code{"jackstraw_rajive"}.
#' @param ... Ignored.
#'
#' @return A \code{data.frame} with columns \code{block}, \code{component},
#'   \code{n_features}, \code{n_significant}, and \code{alpha}.
#'
#' @export
summary.jackstraw_rajive <- function(object, ...) {
  alpha      <- attr(object, "alpha")
  correction <- attr(object, "correction")
  joint_rank <- attr(object, "joint_rank")
  K          <- attr(object, "n_blocks")

  rows <- vector("list", K * joint_rank)
  idx  <- 1L

  for (k in seq_len(K)) {
    for (j in seq_len(joint_rank)) {
      comp   <- object[[k]][[j]]
      n_feat <- length(comp$f_obs)
      n_sig  <- sum(comp$significant, na.rm = TRUE)
      rows[[idx]] <- data.frame(
        block         = names(object)[k],
        component     = paste0("comp", j),
        n_features    = n_feat,
        n_significant = n_sig,
        alpha         = alpha,
        correction    = correction,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  df <- do.call(rbind, rows)
  print(df, row.names = FALSE)
  invisible(df)
}


# ---------------------------------------------------------------------------
# Part 4: Visualization
# ---------------------------------------------------------------------------

#' Plot jackstraw results
#'
#' Produces diagnostic plots for a \code{\link{jackstraw_rajive}} object.
#' Three plot types are supported (selected via \code{type}):
#'
#' \describe{
#'   \item{\code{"pvalue_hist"}}{Histogram of empirical p-values for a single
#'     block / component.  A horizontal reference line is drawn at the
#'     effective significance threshold (Bonferroni-adjusted if applicable).
#'     Enrichment of small p-values indicates real signal.}
#'   \item{\code{"scatter"}}{Scatter plot of observed F-statistic (x-axis)
#'     versus \eqn{-\log_{10}(p\text{-value})} (y-axis) for a single block /
#'     component.  Significant features are coloured red; non-significant
#'     features are grey.  Features are optionally labelled when column names
#'     are present and \code{label_top > 0}.}
#'   \item{\code{"loadings_significance"}}{Heatmap of
#'     \eqn{-\log_{10}(p\text{-value})} across all joint components for a
#'     single block.  Significant cells are marked with an asterisk.}
#' }
#'
#' @param jackstraw_result An object of class \code{"jackstraw_rajive"}.
#' @param type Character string; one of \code{"pvalue_hist"},
#'   \code{"scatter"}, or \code{"loadings_significance"}.
#' @param block Positive integer; which block to plot.  Default \code{1}.
#' @param component Positive integer; which joint component to plot (used for
#'   \code{"pvalue_hist"} and \code{"scatter"}).  Default \code{1}.
#' @param label_top Non-negative integer; for \code{"scatter"}, the number of
#'   top-significant features to label.  Set to \code{0} to suppress labels.
#'   Default \code{10}.
#' @param ... Ignored (reserved for future use).
#'
#' @return A \code{ggplot2} object.
#'
#' @references
#' Yang X, Hoadley KA, Hannig J, Marron JS (2021). Statistical inference for
#' data integration. \emph{arXiv:2109.12272}.
#'
#' @seealso \code{\link{jackstraw_rajive}}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n   <- 50
#' pks <- c(100, 80)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = n,
#'                       pks = pks, dist.type = 1)
#' data.ajive           <- Y$sim_data
#' initial_signal_ranks <- c(5, 4)
#' ajive_result <- Rajive(data.ajive, initial_signal_ranks)
#' js <- jackstraw_rajive(ajive_result, data.ajive, alpha = 0.05, n_null = 10)
#' plot_jackstraw(js, type = "pvalue_hist", block = 1, component = 1)
#' plot_jackstraw(js, type = "scatter",     block = 1, component = 1)
#' plot_jackstraw(js, type = "loadings_significance", block = 1)
#' }
#'
#' @import ggplot2
#' @export
plot_jackstraw <- function(jackstraw_result,
                           type      = c("pvalue_hist", "scatter",
                                         "loadings_significance"),
                           block     = 1L,
                           component = 1L,
                           label_top = 10L,
                           ...) {

  type       <- match.arg(type)
  alpha      <- attr(jackstraw_result, "alpha")
  correction <- attr(jackstraw_result, "correction")
  joint_rank <- attr(jackstraw_result, "joint_rank")
  K          <- attr(jackstraw_result, "n_blocks")

  if (block < 1L || block > K) {
    stop(sprintf("'block' must be between 1 and %d.", K), call. = FALSE)
  }

  block_data <- jackstraw_result[[block]]
  block_nm   <- names(jackstraw_result)[block]

  if (type %in% c("pvalue_hist", "scatter")) {
    if (component < 1L || component > joint_rank) {
      stop(sprintf("'component' must be between 1 and %d.", joint_rank),
           call. = FALSE)
    }

    comp   <- block_data[[component]]
    p_vals <- comp$p_values
    f_obs  <- comp$f_obs
    sig    <- comp$significant
    d      <- length(p_vals)

    # Effective alpha threshold for annotation
    alpha_eff <- if (correction == "bonferroni") {
      alpha / (d * joint_rank)
    } else {
      alpha
    }

    if (type == "pvalue_hist") {
      p_df <- data.frame(p_values = p_vals)

      ggplot(p_df, aes(x = p_values)) +
        geom_histogram(bins = 20, fill = "steelblue", colour = "white",
                       boundary = 0) +
        geom_hline(yintercept = d / 20, linetype = "dashed",
                   colour = "grey50") +
        geom_vline(xintercept = alpha_eff, linetype = "dotted",
                   colour = "red", linewidth = 0.8) +
        labs(
          title    = sprintf("P-value histogram: %s, comp %d", block_nm,
                             component),
          subtitle = sprintf("Red line = significance threshold (%.4g)",
                             alpha_eff),
          x        = "Empirical p-value",
          y        = "Count"
        ) +
        theme_bw()

    } else {
      # "scatter": F-stat vs -log10(p-value)
      feat_names <- names(p_vals)
      if (is.null(feat_names)) feat_names <- as.character(seq_len(d))

      scatter_df <- data.frame(
        feature  = feat_names,
        f_stat   = f_obs,
        neg_logp = -log10(pmax(p_vals, .Machine$double.eps)),
        sig      = ifelse(is.na(sig), FALSE, sig),
        stringsAsFactors = FALSE
      )

      sig_col <- c("FALSE" = "grey70", "TRUE" = "firebrick")

      p <- ggplot(scatter_df,
                  aes(x = f_stat, y = neg_logp,
                      colour = as.character(sig))) +
        geom_point(alpha = 0.7) +
        scale_colour_manual(values = sig_col,
                            labels = c("FALSE" = "Not significant",
                                       "TRUE"  = "Significant"),
                            name   = NULL) +
        geom_hline(yintercept = -log10(alpha_eff), linetype = "dotted",
                   colour = "red", linewidth = 0.8) +
        labs(
          title    = sprintf("F-stat vs -log10(p): %s, comp %d",
                             block_nm, component),
          subtitle = sprintf("Red line = significance threshold (%.4g)",
                             alpha_eff),
          x        = "Observed F-statistic",
          y        = expression(-log[10](p))
        ) +
        theme_bw()

      if (label_top > 0L && any(scatter_df$sig)) {
        top_idx <- order(scatter_df$neg_logp, decreasing = TRUE)[
          seq_len(min(label_top, sum(scatter_df$sig)))
        ]
        lab_df <- scatter_df[top_idx, , drop = FALSE]
        p <- p + ggplot2::geom_text(
          data    = lab_df,
          mapping = aes(label = feature),
          hjust   = -0.1, vjust = 0.5, size = 3,
          colour  = "black"
        )
      }
      p
    }

  } else {
    # "loadings_significance": heatmap over all components for one block
    d <- length(block_data[["comp1"]]$p_values)
    feat_names <- names(block_data[["comp1"]]$p_values)
    if (is.null(feat_names)) feat_names <- as.character(seq_len(d))

    heat_rows <- vector("list", joint_rank)
    for (j in seq_len(joint_rank)) {
      comp   <- block_data[[j]]
      p_vals <- comp$p_values
      sig    <- comp$significant
      d_j    <- length(p_vals)
      alpha_eff <- if (correction == "bonferroni") {
        alpha / (d_j * joint_rank)
      } else {
        alpha
      }

      heat_rows[[j]] <- data.frame(
        feature   = feat_names,
        component = paste0("comp", j),
        neg_logp  = -log10(pmax(p_vals, .Machine$double.eps)),
        sig       = ifelse(is.na(sig), FALSE, sig),
        stringsAsFactors = FALSE
      )
    }

    heat_df <- do.call(rbind, heat_rows)
    heat_df$feature <- factor(heat_df$feature, levels = feat_names)

    p <- ggplot(heat_df,
                aes(x = component, y = feature, fill = neg_logp)) +
      geom_tile(colour = "white") +
      scale_fill_gradient(low = "white", high = "steelblue",
                          name = expression(-log[10](p))) +
      labs(
        title = sprintf("Loadings significance heatmap: %s", block_nm),
        x     = "Joint component",
        y     = "Feature"
      ) +
      theme_bw() +
      theme(axis.text.y = element_text(size = ifelse(d <= 50, 7, 4)))

    # Add asterisks for significant cells
    sig_df <- heat_df[heat_df$sig, , drop = FALSE]
    if (nrow(sig_df) > 0L) {
      p <- p + geom_text(data = sig_df, aes(label = "*"),
                         colour = "black", size = 4, vjust = 0.75)
    }
    p
  }
}


# ---------------------------------------------------------------------------
# Part 5: Convenience accessor
# ---------------------------------------------------------------------------

#' Extract significant variables from jackstraw results
#'
#' Convenience accessor that retrieves the vector of significant variables
#' (either column names or integer indices) for a given block and joint
#' component.
#'
#' @param jackstraw_result An object of class \code{"jackstraw_rajive"}.
#' @param block Positive integer; which block.
#' @param component Positive integer; which joint component.
#'
#' @return A character vector of variable names (when the original data matrix
#'   had column names) or an integer vector of column indices.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n   <- 50
#' pks <- c(100, 80)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = n,
#'                       pks = pks, dist.type = 1)
#' data.ajive           <- Y$sim_data
#' initial_signal_ranks <- c(5, 4)
#' ajive_result <- Rajive(data.ajive, initial_signal_ranks)
#' js <- jackstraw_rajive(ajive_result, data.ajive, alpha = 0.05, n_null = 10)
#' get_significant_vars(js, block = 1, component = 1)
#' }
#'
#' @export
get_significant_vars <- function(jackstraw_result, block = 1L, component = 1L) {
  if (!inherits(jackstraw_result, "jackstraw_rajive")) {
    stop("'jackstraw_result' must be an object of class 'jackstraw_rajive'.",
         call. = FALSE)
  }
  K          <- attr(jackstraw_result, "n_blocks")
  joint_rank <- attr(jackstraw_result, "joint_rank")

  if (block < 1L || block > K) {
    stop(sprintf("'block' must be between 1 and %d.", K), call. = FALSE)
  }
  if (component < 1L || component > joint_rank) {
    stop(sprintf("'component' must be between 1 and %d.", joint_rank),
         call. = FALSE)
  }

  jackstraw_result[[block]][[component]][["significant_vars"]]
}
