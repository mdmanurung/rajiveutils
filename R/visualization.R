# visualization.R — Interpretation and diagnostic helpers for rajiveplus
#
# Public API:
#   extract_components()   — extract tidy data / diagnostic payload from Rajive output
#   rank_features()        — rank features by loading, contribution, overlap, or significance
#   plot_components()      — unified plotting entry point
#   associate_components() — association between component scores and metadata
#   assess_stability()     — bootstrap stability assessment (joint rank / loadings)
#   summarize_components() — summarize outputs across components/diagnostics/tests
#
# Aliases (thin wrappers over rank_features):
#   get_top_loadings()
#   get_feature_contributions()
#   compare_feature_sets_across_blocks()
#   summarize_significant_vars()
#
# Internal helpers (not exported):
#   .extract_rank_diagnostics()
#   .plot_rank_threshold()
#   .plot_bound_distributions()
#   .plot_ajive_diagnostic()
#   .procrustes_align()       — orthogonal Procrustes alignment for loading stability
#   .rf_resolve_method()      — method validation for rank_features
#   .rf_get_loadings()        — loading matrix accessor for rank_features
#   .rf_loading_df()          — top_loadings / contribution worker
#   .rf_overlap()             — overlap worker
#   .rf_significant()         — significant worker


# ---------------------------------------------------------------------------
# extract_components
# ---------------------------------------------------------------------------

#' Extract components or diagnostics from a Rajive decomposition
#'
#' Retrieves structured information from the output of \code{\link{Rajive}}.
#' Currently supports \code{what = "rank_diagnostics"} which returns the joint
#' rank selection diagnostics (observed singular values, Wedin/random-direction
#' bound samples, cutoffs, and identifiability-dropped components).
#'
#' @param ajive_output An object of class \code{"rajive"} (output of
#'   \code{\link{Rajive}}).
#' @param blocks List of original block matrices. Reserved for future
#'   extraction modes; currently unused.
#' @param source Character. One of \code{"rajive"} (default) or
#'   \code{"jackstraw"}.  Selects the result object to extract from.
#'   \code{"jackstraw"} requires \code{jackstraw_result} to be supplied.
#' @param what Character scalar. What to extract.  One of \code{"scores"},
#'   \code{"loadings"}, \code{"variance"}, \code{"significance"}, or
#'   \code{"rank_diagnostics"}.  Only \code{"rank_diagnostics"} is currently
#'   fully implemented.
#' @param type Character. One of \code{"joint"} (default),
#'   \code{"individual"}, or \code{"both"}.  Selects which decomposition
#'   component type to extract.
#' @param block Integer or character. Block identifier (name or 1-based
#'   index) to restrict extraction to a single block.  \code{NULL} (default)
#'   returns results for all blocks.
#' @param component Integer. Component index to extract.  \code{NULL}
#'   (default) returns all components.
#' @param jackstraw_result An optional object of class
#'   \code{"jackstraw_rajive"} (output of \code{\link{jackstraw_rajive}}).
#'   Required when \code{source = "jackstraw"}.
#' @param format Character. \code{"wide"} (default) returns a named list of
#'   class \code{c("rajive_diagnostics", "list")}; \code{"long"} returns a
#'   tidy \code{data.frame} with one row per observed singular value.
#' @param include_meta Logical. If \code{TRUE}, additional metadata fields
#'   useful for plot labels are included in the wide list.
#' @param ... Reserved for future arguments.
#'
#' @return For \code{what = "rank_diagnostics"} and \code{format = "wide"}: a
#'   named list of class \code{c("rajive_diagnostics", "list")} containing:
#'   \describe{
#'     \item{\code{obs_svals}}{Numeric vector of observed common singular values.}
#'     \item{\code{obs_svals_sq}}{Numeric vector of squared common singular values.}
#'     \item{\code{joint_rank_estimate}}{Integer joint rank.}
#'     \item{\code{overall_sv_sq_threshold}}{Numeric combined threshold or
#'       \code{NA_real_} when no bounds were computed.}
#'     \item{\code{wedin_samples}}{Numeric vector of Wedin bound samples, or
#'       \code{NULL}.}
#'     \item{\code{rand_dir_samples}}{Numeric vector of random-direction bound
#'       samples, or \code{NULL}.}
#'     \item{\code{wedin_cutoff}}{5th-percentile Wedin threshold or
#'       \code{NA_real_}.}
#'     \item{\code{rand_cutoff}}{95th-percentile random-direction threshold or
#'       \code{NA_real_}.}
#'     \item{\code{wedin_percentile}}{Percentile used for Wedin cutoff (5) or
#'       \code{NA_real_}.}
#'     \item{\code{rand_percentile}}{Percentile used for random cutoff (95) or
#'       \code{NA_real_}.}
#'     \item{\code{identif_dropped}}{Integer vector of component indices
#'       dropped by the identifiability check (length 0 when none dropped).}
#'     \item{\code{cutoff_rule}}{One of \code{"max(wedin, random)"},
#'       \code{"wedin_only"}, \code{"random_only"}, \code{"none_available"}.}
#'     \item{\code{has_wedin}}{Logical.}
#'     \item{\code{has_random}}{Logical.}
#'   }
#'   For \code{format = "long"}: a \code{data.frame} with columns
#'   \code{component_index}, \code{obs_sval}, \code{obs_sval_sq},
#'   \code{classification}, \code{joint_rank_estimate},
#'   \code{overall_sv_sq_threshold}, \code{wedin_cutoff}, \code{rand_cutoff}.
#'
#' @examples
#' \donttest{
#' n   <- 40; pks <- c(30, 20)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#' diag_wide <- extract_components(fit, what = "rank_diagnostics")
#' diag_long <- extract_components(fit, what = "rank_diagnostics", format = "long")
#' }
#'
#' @export
extract_components <- function(ajive_output = NULL,
                               blocks = NULL,
                               source = c("rajive", "jackstraw"),
                               what   = c("scores", "loadings", "variance",
                                          "significance", "rank_diagnostics"),
                               type = c("joint", "individual", "both"),
                               block = NULL,
                               component = NULL,
                               jackstraw_result = NULL,
                               format = c("wide", "long"),
                               include_meta = FALSE,
                               ...) {
  source <- match.arg(source)
  what   <- match.arg(what)
  type   <- match.arg(type)
  format <- match.arg(format)

  if (what == "rank_diagnostics") {
    if (!inherits(ajive_output, "rajive")) {
      stop("ajive_output must be of class 'rajive'.", call. = FALSE)
    }
    return(.extract_rank_diagnostics(ajive_output, format = format,
                                     include_meta = include_meta))
  }

  if (what %in% c("scores", "loadings", "variance")) {
    if (!inherits(ajive_output, "rajive")) {
      stop("ajive_output must be of class 'rajive'.", call. = FALSE)
    }
  }

  if (what == "scores") {
    return(.extract_scores(ajive_output, type = type, block = block,
                           component = component, format = format))
  }

  if (what == "loadings") {
    return(.extract_loadings(ajive_output, type = type, block = block,
                             component = component, format = format))
  }

  if (what == "variance") {
    if (is.null(blocks) || !is.list(blocks)) {
      stop("`blocks` must be supplied as a list for what = 'variance'.",
           call. = FALSE)
    }
    out <- as.data.frame(showVarExplained_robust(ajive_output, blocks),
                         stringsAsFactors = FALSE)
    if (!is.null(block) && "Block" %in% names(out)) {
      b <- as.integer(block)
      out <- out[out$Block %in% b | out$Block %in% paste0("block", b), , drop = FALSE]
    }
    return(out)
  }

  if (what == "significance") {
    if (!identical(source, "jackstraw") && is.null(jackstraw_result)) {
      stop("what = 'significance' requires source = 'jackstraw' or jackstraw_result.",
           call. = FALSE)
    }
    if (is.null(jackstraw_result)) {
      stop("`jackstraw_result` must be supplied for what = 'significance'.",
           call. = FALSE)
    }
    return(rank_features(jackstraw_result = jackstraw_result,
                         mode = "significant",
                         block = block,
                         component = component))
  }

  stop("Unsupported extraction request.", call. = FALSE)
}


# ---------------------------------------------------------------------------
# plot_components
# ---------------------------------------------------------------------------

#' Plot AJIVE diagnostic visualizations
#'
#' Produces diagnostic plots from the output of \code{\link{Rajive}}.  Three
#' plot types are supported:
#' \describe{
#'   \item{\code{"rank_threshold"}}{Squared observed singular values with
#'     Wedin and/or random-direction cutoff lines, colored by joint/nonjoint
#'     classification.  The active threshold is the
#'     \code{max(wedin, random_direction)} rule -- a practical heuristic that
#'     is conservative in practice but does not carry formal FWER or FDR
#'     control for rank selection (see \code{StatisticalAudits.md}, Finding 6).}
#'   \item{\code{"bound_distributions"}}{Histogram(s) of Wedin and/or
#'     random-direction bound samples with percentile cutoff lines.}
#'   \item{\code{"ajive_diagnostic"}}{Full composite panel combining the
#'     threshold plot (top) with bound-distribution histograms (bottom).
#'     Requires the \pkg{patchwork} package.}
#' }
#'
#' @param ajive_output An object of class \code{"rajive"} (output of
#'   \code{\link{Rajive}}), or \code{NULL} when a precomputed diagnostics
#'   payload is supplied via \code{diagnostics}.
#' @param blocks List of original block matrices. Used by score/loading plot
#'   types that require the raw data; ignored for diagnostic plots.
#' @param jackstraw_result An optional object of class
#'   \code{"jackstraw_rajive"} (output of \code{\link{jackstraw_rajive}}).
#'   When supplied, significance overlays are added to supported plot types.
#' @param assoc_results Optional output of \code{\link{associate_components}}.
#'   When supplied, association plots use precomputed results instead of
#'   re-computing them.
#' @param stability_result Optional output of \code{\link{assess_stability}}.
#'   Required for \code{plot_type = "stability"}.
#' @param plot_type Character.  One of \code{"ajive_diagnostic"},
#'   \code{"rank_threshold"}, \code{"bound_distributions"},
#'   \code{"pairs"}, \code{"density"}, \code{"top_features"},
#'   \code{"component_heatmap"}, \code{"variance"}, \code{"association"},
#'   \code{"volcano"}, \code{"jackstraw_summary"}, or \code{"stability"}.
#' @param type Character. One of \code{"joint"} (default),
#'   \code{"individual"}, or \code{"both"}.  Selects which decomposition
#'   component type to plot.
#' @param block Integer or character. Block identifier (name or 1-based
#'   index) to restrict the plot to a single data block.  \code{NULL}
#'   (default) plots all blocks where applicable.
#' @param component Integer. Component index to plot.  \code{NULL} (default)
#'   plots all components where applicable.
#' @param group Optional grouping vector (length equal to the number of
#'   samples) or a column name in metadata used to colour points.
#' @param style Optional named list of styling overrides (e.g.,
#'   \code{list(point_size = 1.5, alpha = 0.6)}) passed through to the
#'   underlying plot helpers.
#' @param top_n Integer. Number of top features to display in feature-ranking
#'   plot types (default \code{20L}).
#' @param ... Additional arguments passed to the plot helpers:
#'   \describe{
#'     \item{\code{diagnostics}}{Precomputed diagnostics list from
#'       \code{extract_components(..., what = "rank_diagnostics")}.  When
#'       supplied, \code{ajive_output} is ignored.}
#'     \item{\code{theme_fn}}{A \code{ggplot2} theme function; defaults to
#'       \code{ggplot2::theme_bw}.}
#'     \item{\code{cutoff_quantile}}{Numeric in (0, 1). Override the
#'       percentile cutoff for bound-distribution plots (default: use the
#'       cutoffs stored in the diagnostics payload).}
#'   }
#'
#' @return A \code{ggplot} object (or a \pkg{patchwork} composite for
#'   \code{plot_type = "ajive_diagnostic"}).
#'
#' @section Variance explained language:
#'   When variance-explained plots are added in a future phase, outputs will
#'   distinguish \emph{joint} variance explained from
#'   \emph{block-specific individual} variance explained after joint signal
#'   removal.  The two quantities are not interchangeable and should not be
#'   described with generic PCA-style language (StatisticalAudits.md, Finding 7).
#'
#' @examples
#' \donttest{
#' n   <- 40; pks <- c(30, 20)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#' plot_components(fit, plot_type = "rank_threshold")
#' plot_components(fit, plot_type = "bound_distributions")
#' plot_components(fit, plot_type = "ajive_diagnostic")
#' }
#'
#' @export
plot_components <- function(ajive_output = NULL,
                            blocks = NULL,
                            jackstraw_result = NULL,
                            assoc_results = NULL,
                            stability_result = NULL,
                            plot_type = c("pairs", "density", "top_features",
                                          "component_heatmap", "variance",
                                          "association", "volcano",
                                          "jackstraw_summary", "stability",
                                          "ajive_diagnostic", "rank_threshold",
                                          "bound_distributions"),
                            type = c("joint", "individual", "both"),
                            block = NULL,
                            component = NULL,
                            group = NULL,
                            style = NULL,
                            top_n = 20L,
                            ...) {
  plot_type <- match.arg(plot_type)
  type <- match.arg(type)
  dots <- list(...)

  if (plot_type %in% c("ajive_diagnostic", "rank_threshold", "bound_distributions")) {
    diag <- if (!is.null(dots$diagnostics)) {
      dots$diagnostics
    } else if (!is.null(ajive_output)) {
      if (!inherits(ajive_output, "rajive")) {
        stop("ajive_output must be of class 'rajive'.", call. = FALSE)
      }
      extract_components(ajive_output, what = "rank_diagnostics", format = "wide")
    } else {
      stop("ajive_output must be provided when diagnostics is not supplied.",
           call. = FALSE)
    }

    return(switch(plot_type,
      ajive_diagnostic    = .plot_ajive_diagnostic(diag, dots),
      rank_threshold      = .plot_rank_threshold(diag, dots),
      bound_distributions = .plot_bound_distributions(diag, dots)
    ))
  }

  if (plot_type == "pairs") {
    return(.plot_pairs_scores(ajive_output, type = type, block = block,
                              group = group, component = component))
  }
  if (plot_type == "density") {
    return(.plot_density_scores(ajive_output, type = type, block = block,
                                group = group, component = component))
  }
  if (plot_type == "top_features") {
    return(.plot_top_features(ajive_output, type = type, block = block,
                              component = component, top_n = as.integer(top_n)))
  }
  if (plot_type == "component_heatmap") {
    return(.plot_component_heatmap(ajive_output, type = type, block = block,
                                   component = component))
  }
  if (plot_type == "variance") {
    if (is.null(blocks))
      stop("`blocks` must be supplied for plot_type = 'variance'.", call. = FALSE)
    return(.plot_variance_summary(ajive_output, blocks = blocks))
  }
  if (plot_type == "association") {
    return(.plot_association_summary(assoc_results = assoc_results, style = style))
  }
  if (plot_type == "volcano") {
    return(.plot_jackstraw_volcano(jackstraw_result = jackstraw_result,
                                   block = block, component = component))
  }
  if (plot_type == "jackstraw_summary") {
    return(.plot_jackstraw_summary(jackstraw_result = jackstraw_result))
  }
  if (plot_type == "stability") {
    return(.plot_stability_summary(stability_result = stability_result, style = style))
  }

  stop("Unsupported plot_type.", call. = FALSE)
}


# ---------------------------------------------------------------------------
# Internal: .extract_rank_diagnostics
# ---------------------------------------------------------------------------

.extract_rank_diagnostics <- function(ajive_output, format, include_meta) {
  jrs <- ajive_output$joint_rank_sel
  if (is.null(jrs)) {
    stop(paste0(
      "joint_rank_sel not found in ajive_output. ",
      "Ensure the object was created by Rajive() without suppressing rank selection."),
      call. = FALSE)
  }

  obs_svals <- jrs[["obs_svals"]]
  if (is.null(obs_svals) || length(obs_svals) == 0L) {
    stop("Observed singular values not found in joint_rank_sel.", call. = FALSE)
  }
  obs_svals_sq <- obs_svals^2

  joint_rank_estimate <- as.integer(ajive_output$joint_rank)

  # Wedin bound -----------------------------------------------------------
  has_wedin <- !is.null(jrs[["wedin"]])
  if (has_wedin) {
    wedin_cutoff     <- as.numeric(jrs[["wedin"]][["wedin_svsq_threshold"]])
    wedin_samples    <- as.numeric(jrs[["wedin"]][["wedin_samples"]])
    wedin_percentile <- 5  # percentile used by get_joint_scores_robustH
  } else {
    wedin_cutoff     <- NA_real_
    wedin_samples    <- NULL
    wedin_percentile <- NA_real_
  }

  # Random direction bound ------------------------------------------------
  has_random <- !is.null(jrs[["rand_dir"]])
  if (has_random) {
    rand_cutoff      <- as.numeric(jrs[["rand_dir"]][["rand_dir_svsq_threshold"]])
    rand_dir_samples <- as.numeric(jrs[["rand_dir"]][["rand_dir_samples"]])
    rand_percentile  <- 95  # percentile used by get_joint_scores_robustH
  } else {
    rand_cutoff      <- NA_real_
    rand_dir_samples <- NULL
    rand_percentile  <- NA_real_
  }

  # Permutation bound -----------------------------------------------------
  has_perm <- !is.null(jrs[["perm"]])
  if (has_perm) {
    perm_cutoff   <- as.numeric(jrs[["perm"]][["perm_svsq_threshold"]])
    perm_samples  <- as.numeric(jrs[["perm"]][["perm_samples"]])
    perm_percentile <- 95
  } else {
    perm_cutoff  <- NA_real_
    perm_samples <- NULL
    perm_percentile <- NA_real_
  }

  if (!has_wedin && !has_random && !has_perm) {
    warning(paste0(
      "Wedin, random-direction, and permutation bounds are absent from joint_rank_sel. ",
      "Diagnostic visualizations will show observed singular values only."),
      call. = FALSE)
  }

  overall_sv_sq_threshold <- jrs[["overall_sv_sq_threshold"]]
  if (is.null(overall_sv_sq_threshold)) overall_sv_sq_threshold <- NA_real_

  cutoff_rule <- if (has_wedin && has_random) {
    "max(wedin, random)"
  } else if (has_wedin && has_perm) {
    "max(wedin, perm)"
  } else if (has_wedin) {
    "wedin_only"
  } else if (has_perm) {
    "perm_only"
  } else if (has_random) {
    "random_only"
  } else {
    "none_available"
  }

  identif_dropped <- jrs[["identif_dropped"]]
  if (is.null(identif_dropped)) identif_dropped <- integer(0L)

  payload <- list(
    obs_svals               = obs_svals,
    obs_svals_sq            = obs_svals_sq,
    joint_rank_estimate     = joint_rank_estimate,
    overall_sv_sq_threshold = overall_sv_sq_threshold,
    wedin_samples           = wedin_samples,
    rand_dir_samples        = rand_dir_samples,
    perm_samples            = perm_samples,
    wedin_cutoff            = wedin_cutoff,
    rand_cutoff             = rand_cutoff,
    perm_cutoff             = perm_cutoff,
    wedin_percentile        = wedin_percentile,
    rand_percentile         = rand_percentile,
    perm_percentile         = perm_percentile,
    identif_dropped         = identif_dropped,
    cutoff_rule             = cutoff_rule,
    has_wedin               = has_wedin,
    has_random              = has_random,
    has_perm                = has_perm
  )
  class(payload) <- c("rajive_diagnostics", "list")

  if (format == "wide") return(payload)

  # Long format: one row per observed singular value ----------------------
  n_svals <- length(obs_svals)
  is_joint <- seq_len(n_svals) <= joint_rank_estimate
  classification <- ifelse(
    seq_len(n_svals) %in% identif_dropped, "dropped",
    ifelse(is_joint, "joint", "nonjoint")
  )

  data.frame(
    component_index         = seq_len(n_svals),
    obs_sval                = obs_svals,
    obs_sval_sq             = obs_svals_sq,
    classification          = classification,
    joint_rank_estimate     = joint_rank_estimate,
    overall_sv_sq_threshold = overall_sv_sq_threshold,
    wedin_cutoff            = wedin_cutoff,
    rand_cutoff             = rand_cutoff,
    perm_cutoff             = perm_cutoff,
    stringsAsFactors        = FALSE
  )
}


# Internal: score/loading extractors for extract_components
.extract_scores <- function(ajive_output, type, block, component, format) {
  K <- length(ajive_output$block_decomps) %/% 3L
  blocks_k <- if (is.null(block)) seq_len(K) else as.integer(block)
  use_types <- if (type == "both") c("joint", "individual") else type

  out_list <- list()
  long_rows <- list()

  for (k in blocks_k) {
    if (k < 1L || k > K) next
    for (tp in use_types) {
      mat <- if (tp == "joint") {
        ajive_output$joint_scores
      } else {
        ajive_output$block_decomps[[3L * (k - 1L) + 1L]]$u
      }
      if (is.null(mat) || ncol(mat) == 0L) next
      comps <- if (is.null(component)) seq_len(ncol(mat)) else as.integer(component)
      comps <- comps[comps >= 1L & comps <= ncol(mat)]
      if (length(comps) == 0L) next
      mat <- mat[, comps, drop = FALSE]

      key <- paste0("block", k, "_", tp, "_scores")
      out_list[[key]] <- mat

      if (format == "long") {
        df <- as.data.frame(mat, stringsAsFactors = FALSE)
        colnames(df) <- paste0("comp", comps)
        df$sample <- seq_len(nrow(df))
        df_long <- stats::reshape(
          df,
          varying = paste0("comp", comps),
          v.names = "score",
          timevar = "component",
          times = comps,
          direction = "long"
        )
        rownames(df_long) <- NULL
        df_long$component <- as.integer(sub("^comp", "", as.character(df_long$component)))
        df_long$block <- k
        df_long$type <- tp
        long_rows[[length(long_rows) + 1L]] <-
          df_long[, c("block", "type", "component", "sample", "score")]
      }
    }
  }

  if (format == "wide") return(out_list)
  if (length(long_rows) == 0L) {
    return(data.frame(block = integer(0), type = character(0), component = integer(0),
                      sample = integer(0), score = numeric(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, long_rows)
}

.extract_loadings <- function(ajive_output, type, block, component, format) {
  K <- length(ajive_output$block_decomps) %/% 3L
  blocks_k <- if (is.null(block)) seq_len(K) else as.integer(block)
  use_types <- if (type == "both") c("joint", "individual") else type

  out_list <- list()
  long_rows <- list()

  for (k in blocks_k) {
    if (k < 1L || k > K) next
    for (tp in use_types) {
      idx <- if (tp == "individual") 3L * (k - 1L) + 1L else 3L * (k - 1L) + 2L
      mat <- ajive_output$block_decomps[[idx]]$v
      if (is.null(mat) || ncol(mat) == 0L) next
      comps <- if (is.null(component)) seq_len(ncol(mat)) else as.integer(component)
      comps <- comps[comps >= 1L & comps <= ncol(mat)]
      if (length(comps) == 0L) next
      mat <- mat[, comps, drop = FALSE]

      key <- paste0("block", k, "_", tp, "_loadings")
      out_list[[key]] <- mat

      if (format == "long") {
        df <- as.data.frame(mat, stringsAsFactors = FALSE)
        colnames(df) <- paste0("comp", comps)
        df$feature <- seq_len(nrow(df))
        df_long <- stats::reshape(
          df,
          varying = paste0("comp", comps),
          v.names = "loading",
          timevar = "component",
          times = comps,
          direction = "long"
        )
        rownames(df_long) <- NULL
        df_long$component <- as.integer(sub("^comp", "", as.character(df_long$component)))
        df_long$block <- k
        df_long$type <- tp
        long_rows[[length(long_rows) + 1L]] <-
          df_long[, c("block", "type", "component", "feature", "loading")]
      }
    }
  }

  if (format == "wide") return(out_list)
  if (length(long_rows) == 0L) {
    return(data.frame(block = integer(0), type = character(0), component = integer(0),
                      feature = integer(0), loading = numeric(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, long_rows)
}


# ---------------------------------------------------------------------------
# Internal: .plot_rank_threshold
# ---------------------------------------------------------------------------

.plot_rank_threshold <- function(diag, dots) {
  theme_fn <- if (!is.null(dots$theme_fn)) dots$theme_fn else ggplot2::theme_bw

  n_svals <- length(diag$obs_svals)
  classification <- ifelse(
    seq_len(n_svals) %in% diag$identif_dropped, "dropped",
    ifelse(seq_len(n_svals) <= diag$joint_rank_estimate, "joint", "nonjoint")
  )

  df <- data.frame(
    idx     = seq_len(n_svals),
    sval_sq = diag$obs_svals_sq,
    class   = factor(classification, levels = c("joint", "nonjoint", "dropped"))
  )

  colors <- c(joint = "black",  nonjoint = "grey60", dropped = "grey30")
  shapes <- c(joint = 16L, nonjoint = 1L,  dropped = 4L)

  has_wedin  <- diag$has_wedin
  has_random <- diag$has_random
  has_perm   <- isTRUE(diag$has_perm)
  threshold  <- diag$overall_sv_sq_threshold

  subtitle <- if (has_wedin && has_random) {
    sprintf("Joint rank: %d | Threshold: %.4f (%s)",
            diag$joint_rank_estimate, threshold, diag$cutoff_rule)
  } else if (has_wedin && has_perm) {
    sprintf("Joint rank: %d | Threshold: %.4f (%s)",
            diag$joint_rank_estimate, threshold, diag$cutoff_rule)
  } else if (has_wedin) {
    sprintf("Joint rank: %d | Wedin bound only (rand_dir/perm bound unavailable)",
            diag$joint_rank_estimate)
  } else if (has_perm) {
    sprintf("Joint rank: %d | Permutation bound only (Wedin bound unavailable)",
            diag$joint_rank_estimate)
  } else if (has_random) {
    sprintf("Joint rank: %d | Random direction bound only (Wedin bound unavailable)",
            diag$joint_rank_estimate)
  } else {
    sprintf("Joint rank: %d | diagnostic bounds unavailable",
            diag$joint_rank_estimate)
  }

  p <- ggplot2::ggplot(df,
      ggplot2::aes(x = .data$idx, y = .data$sval_sq,
                   color = .data$class, shape = .data$class)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = colors, name = "Classification",
                                drop = FALSE) +
    ggplot2::scale_shape_manual(values = shapes, name = "Classification",
                                drop = FALSE) +
    ggplot2::labs(
      x        = "Component index",
      y        = expression(Squared~singular~value~(sigma^2)),
      title    = "AJIVE Joint Rank Threshold",
      subtitle = subtitle
    ) +
    theme_fn()

  if (has_wedin) {
    wedin_lw <- if (!has_random || diag$wedin_cutoff >= diag$rand_cutoff) 1.5 else 0.75
    p <- p + ggplot2::geom_hline(
      yintercept = diag$wedin_cutoff,
      color      = "steelblue",
      linetype   = "dashed",
      linewidth  = wedin_lw
    )
  }

  if (has_random) {
    rand_lw <- if (!has_wedin || diag$rand_cutoff >= diag$wedin_cutoff) 1.5 else 0.75
    p <- p + ggplot2::geom_hline(
      yintercept = diag$rand_cutoff,
      color      = "tomato",
      linetype   = "dotted",
      linewidth  = rand_lw
    )
  }

  if (has_perm) {
    perm_lw <- if (!has_wedin || diag$perm_cutoff >= diag$wedin_cutoff) 1.5 else 0.75
    p <- p + ggplot2::geom_hline(
      yintercept = diag$perm_cutoff,
      color      = "forestgreen",
      linetype   = "dotdash",
      linewidth  = perm_lw
    )
  }

  p
}


# ---------------------------------------------------------------------------
# Internal: .plot_bound_distributions
# ---------------------------------------------------------------------------

.plot_bound_distributions <- function(diag, dots) {
  theme_fn <- if (!is.null(dots$theme_fn)) dots$theme_fn else ggplot2::theme_bw

  if (!is.null(dots$cutoff_quantile)) {
    q <- dots$cutoff_quantile
    if (!is.numeric(q) || length(q) != 1L || q <= 0 || q >= 1) {
      stop("cutoff_quantile must be in (0, 1).", call. = FALSE)
    }
  }

  plots <- list()

  if (diag$has_wedin && !is.null(diag$wedin_samples)) {
    wedin_cutoff <- if (!is.null(dots$cutoff_quantile)) {
      stats::quantile(diag$wedin_samples, dots$cutoff_quantile)
    } else {
      diag$wedin_cutoff
    }
    df_w <- data.frame(sample = diag$wedin_samples)
    plots[["wedin"]] <- ggplot2::ggplot(df_w, ggplot2::aes(x = .data$sample)) +
      ggplot2::geom_histogram(bins = 30L, fill = "steelblue", alpha = 0.7,
                              color = "white") +
      ggplot2::geom_vline(xintercept = wedin_cutoff, color = "steelblue",
                          linewidth = 1.2, linetype = "dashed") +
      ggplot2::labs(
        x     = "Wedin bound samples",
        y     = "Count",
        title = sprintf("Wedin bound (cutoff = %.4f)", wedin_cutoff)
      ) +
      theme_fn()
  }

  if (diag$has_random && !is.null(diag$rand_dir_samples)) {
    rand_cutoff <- if (!is.null(dots$cutoff_quantile)) {
      stats::quantile(diag$rand_dir_samples, dots$cutoff_quantile)
    } else {
      diag$rand_cutoff
    }
    df_r <- data.frame(sample = diag$rand_dir_samples)
    plots[["rand_dir"]] <- ggplot2::ggplot(df_r, ggplot2::aes(x = .data$sample)) +
      ggplot2::geom_histogram(bins = 30L, fill = "tomato", alpha = 0.7,
                              color = "white") +
      ggplot2::geom_vline(xintercept = rand_cutoff, color = "tomato",
                          linewidth = 1.2, linetype = "dashed") +
      ggplot2::labs(
        x     = "Random direction bound samples",
        y     = "Count",
        title = sprintf("Random direction bound (cutoff = %.4f)", rand_cutoff)
      ) +
      theme_fn()
  }

  if (isTRUE(diag$has_perm) && !is.null(diag$perm_samples)) {
    perm_cutoff <- if (!is.null(dots$cutoff_quantile)) {
      stats::quantile(diag$perm_samples, dots$cutoff_quantile)
    } else {
      diag$perm_cutoff
    }
    df_p <- data.frame(sample = diag$perm_samples)
    plots[["perm"]] <- ggplot2::ggplot(df_p, ggplot2::aes(x = .data$sample)) +
      ggplot2::geom_histogram(bins = 30L, fill = "forestgreen", alpha = 0.7,
                              color = "white") +
      ggplot2::geom_vline(xintercept = perm_cutoff, color = "forestgreen",
                          linewidth = 1.2, linetype = "dashed") +
      ggplot2::labs(
        x     = "Permutation bound samples",
        y     = "Count",
        title = sprintf("Permutation bound (cutoff = %.4f)", perm_cutoff)
      ) +
      theme_fn()
  }

  if (length(plots) == 0L) {
    # No distributions available — annotated empty plot
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label  = "diagnostic bounds unavailable",
                          size   = 5, color = "grey50") +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::theme_void()
    )
  }

  if (length(plots) == 1L) return(plots[[1L]])

  # Multiple distributions: combine side-by-side
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(paste0(
      "Package 'patchwork' is required to combine bound-distribution plots. ",
      "Install with: install.packages('patchwork')."),
      call. = FALSE)
  }
  Reduce(`+`, plots)
}


# ---------------------------------------------------------------------------
# Internal: .plot_ajive_diagnostic
# ---------------------------------------------------------------------------

.plot_ajive_diagnostic <- function(diag, dots) {
  p_thresh <- .plot_rank_threshold(diag, dots)

  has_wedin  <- diag$has_wedin
  has_random <- diag$has_random
  has_perm   <- isTRUE(diag$has_perm)

  # No bounds at all: return single panel with annotation already in subtitle
  if (!has_wedin && !has_random && !has_perm) {
    return(p_thresh)
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(paste0(
      "Package 'patchwork' is required for the composite diagnostic panel. ",
      "Install with: install.packages('patchwork')."),
      call. = FALSE)
  }

  # Build individual distribution plots so we can arrange them
  theme_fn <- if (!is.null(dots$theme_fn)) dots$theme_fn else ggplot2::theme_bw

  resolve_cutoff <- function(samples, stored_cutoff, q) {
    if (!is.null(q)) stats::quantile(samples, q) else stored_cutoff
  }
  q <- dots$cutoff_quantile

  dist_plots <- list()

  if (has_wedin && !is.null(diag$wedin_samples)) {
    wc <- resolve_cutoff(diag$wedin_samples, diag$wedin_cutoff, q)
    df_w <- data.frame(sample = diag$wedin_samples)
    dist_plots[["wedin"]] <- ggplot2::ggplot(df_w, ggplot2::aes(x = .data$sample)) +
      ggplot2::geom_histogram(bins = 30L, fill = "steelblue", alpha = 0.7,
                              color = "white") +
      ggplot2::geom_vline(xintercept = wc, color = "steelblue",
                          linewidth = 1.2, linetype = "dashed") +
      ggplot2::labs(x = "Wedin samples", y = "Count",
                    title = sprintf("Wedin (cutoff = %.4f)", wc)) +
      theme_fn()
  }

  if (has_random && !is.null(diag$rand_dir_samples)) {
    rc <- resolve_cutoff(diag$rand_dir_samples, diag$rand_cutoff, q)
    df_r <- data.frame(sample = diag$rand_dir_samples)
    dist_plots[["rand_dir"]] <- ggplot2::ggplot(df_r, ggplot2::aes(x = .data$sample)) +
      ggplot2::geom_histogram(bins = 30L, fill = "tomato", alpha = 0.7,
                              color = "white") +
      ggplot2::geom_vline(xintercept = rc, color = "tomato",
                          linewidth = 1.2, linetype = "dashed") +
      ggplot2::labs(x = "Rand. dir. samples", y = "Count",
                    title = sprintf("Rand. dir. (cutoff = %.4f)", rc)) +
      theme_fn()
  }

  if (has_perm && !is.null(diag$perm_samples)) {
    pc <- resolve_cutoff(diag$perm_samples, diag$perm_cutoff, q)
    df_p <- data.frame(sample = diag$perm_samples)
    dist_plots[["perm"]] <- ggplot2::ggplot(df_p, ggplot2::aes(x = .data$sample)) +
      ggplot2::geom_histogram(bins = 30L, fill = "forestgreen", alpha = 0.7,
                              color = "white") +
      ggplot2::geom_vline(xintercept = pc, color = "forestgreen",
                          linewidth = 1.2, linetype = "dashed") +
      ggplot2::labs(x = "Permutation samples", y = "Count",
                    title = sprintf("Permutation (cutoff = %.4f)", pc)) +
      theme_fn()
  }

  n_dist <- length(dist_plots)
  bottom <- if (n_dist == 0L) {
    NULL
  } else if (n_dist == 1L) {
    dist_plots[[1L]]
  } else {
    Reduce(`+`, dist_plots)
  }

  if (is.null(bottom)) return(p_thresh)
  p_thresh / bottom
}


# ---------------------------------------------------------------------------
# Additional plot_components workers
# ---------------------------------------------------------------------------

.resolve_plot_scores <- function(ajive_output, type, block) {
  if (!inherits(ajive_output, "rajive"))
    stop("`ajive_output` must be an object of class 'rajive'.", call. = FALSE)

  if (type == "joint") {
    sc <- ajive_output$joint_scores
    if (is.null(sc)) stop("Joint scores are unavailable.", call. = FALSE)
    return(sc)
  }

  if (is.null(block))
    stop("`block` must be supplied for individual score plotting.", call. = FALSE)
  k <- as.integer(block)[1L]
  idx <- 3L * (k - 1L) + 1L
  sc <- ajive_output$block_decomps[[idx]]$u
  if (is.null(sc)) stop("Individual scores are unavailable for requested block.", call. = FALSE)
  sc
}

.plot_pairs_scores <- function(ajive_output, type, block, group, component) {
  sc <- .resolve_plot_scores(ajive_output, if (type == "both") "joint" else type, block)
  comps <- if (is.null(component)) seq_len(min(2L, ncol(sc))) else as.integer(component)
  if (length(comps) < 2L) {
    # Only 1 component available — fall back to component vs sample-index scatter
    message("Only 1 component available; plotting component 1 vs sample index.")
    df2 <- data.frame(x = seq_len(nrow(sc)), y = sc[, 1L], stringsAsFactors = FALSE)
    if (!is.null(group) && length(group) == nrow(df2)) df2$group <- as.factor(group)
    aes2 <- if ("group" %in% names(df2)) ggplot2::aes(x = .data$x, y = .data$y, color = .data$group) else ggplot2::aes(x = .data$x, y = .data$y)
    return(ggplot2::ggplot(df2, aes2) +
             ggplot2::geom_point(alpha = 0.8) +
             ggplot2::theme_bw() +
             ggplot2::labs(x = "Sample index", y = "Component 1 score",
                           title = "Component 1 scores (only 1 component available)"))
  }
  comps <- comps[1:2]
  df <- data.frame(x = sc[, comps[1L]], y = sc[, comps[2L]], stringsAsFactors = FALSE)
  if (!is.null(group) && length(group) == nrow(df)) df$group <- as.factor(group)
  aes_map <- if ("group" %in% names(df)) {
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$group)
  } else {
    ggplot2::aes(x = .data$x, y = .data$y)
  }
  ggplot2::ggplot(df, aes_map) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = paste0("Component ", comps[1L]),
                  y = paste0("Component ", comps[2L]),
                  title = "Component score pairs")
}

.plot_density_scores <- function(ajive_output, type, block, group, component) {
  sc <- .resolve_plot_scores(ajive_output, if (type == "both") "joint" else type, block)
  comp <- if (is.null(component)) 1L else as.integer(component)[1L]
  if (comp < 1L || comp > ncol(sc)) stop("Requested component is out of range.", call. = FALSE)
  df <- data.frame(score = sc[, comp], stringsAsFactors = FALSE)
  if (!is.null(group) && length(group) == nrow(df)) df$group <- as.factor(group)
  if ("group" %in% names(df)) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data$score, fill = .data$group)) +
      ggplot2::geom_density(alpha = 0.35) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = paste0("Component ", comp, " score"),
                    y = "Density", title = "Score density by group")
  } else {
    ggplot2::ggplot(df, ggplot2::aes(x = .data$score)) +
      ggplot2::geom_density(fill = "grey70", alpha = 0.6) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = paste0("Component ", comp, " score"),
                    y = "Density", title = "Score density")
  }
}

.plot_top_features <- function(ajive_output, type, block, component, top_n) {
  tp <- if (type == "both") "joint" else type
  rf <- rank_features(ajive_output = ajive_output,
                      mode = "top_loadings",
                      type = tp,
                      block = block,
                      component = component,
                      n = top_n)
  if (is.null(rf) || nrow(rf) == 0L)
    stop("No ranked features available for plotting.", call. = FALSE)
  rf$feature_label <- ifelse(is.na(rf$feature_name),
                             paste0("f", rf$feature_index), rf$feature_name)
  rf <- rf[order(rf$rank), , drop = FALSE]
  ggplot2::ggplot(rf,
                  ggplot2::aes(x = stats::reorder(.data$feature_label, .data$abs_loading),
                               y = .data$abs_loading)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Feature", y = "Absolute loading", title = "Top-ranked features")
}

.plot_component_heatmap <- function(ajive_output, type, block, component) {
  sc <- .resolve_plot_scores(ajive_output, if (type == "both") "joint" else type, block)
  comps <- if (is.null(component)) seq_len(ncol(sc)) else as.integer(component)
  comps <- comps[comps >= 1L & comps <= ncol(sc)]
  if (length(comps) == 0L) stop("No valid components selected.", call. = FALSE)
  mat <- sc[, comps, drop = FALSE]
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(df) <- c("sample", "component", "value")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$component, y = .data$sample, fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Component", y = "Sample", title = "Component score heatmap")
}

.plot_variance_summary <- function(ajive_output, blocks) {
  vt <- as.data.frame(showVarExplained_robust(ajive_output, blocks), stringsAsFactors = FALSE)
  nms <- names(vt)
  block_col <- nms[1L]
  long <- stats::reshape(vt,
                        varying = setdiff(names(vt), block_col),
                        v.names = "percent",
                        timevar = "part",
                        times = setdiff(names(vt), block_col),
                        direction = "long")
  rownames(long) <- NULL
  names(long)[1L] <- "block"
  ggplot2::ggplot(long, ggplot2::aes(x = .data$block, y = .data$percent, fill = .data$part)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Block", y = "Percent", fill = "Part", title = "Variance explained")
}

.plot_association_summary <- function(assoc_results, style = NULL) {
  if (is.null(assoc_results) || !is.data.frame(assoc_results))
    stop("`assoc_results` must be provided for plot_type = 'association'.", call. = FALSE)
  req <- c("variable", "component", "p_adj")
  miss <- setdiff(req, names(assoc_results))
  if (length(miss) > 0L) stop("assoc_results missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

  st <- if (is.null(style)) "heatmap" else style
  if (identical(st, "forest")) {
    d <- assoc_results
    d$neglog10 <- -log10(pmax(d$p_adj, .Machine$double.eps))
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = .data$neglog10,
                                      y = paste0(.data$variable, " (C", .data$component, ")"))) +
        ggplot2::geom_point() + ggplot2::theme_bw() +
        ggplot2::labs(x = "-log10(adj p)", y = "Association", title = "Association forest")
    )
  }

  d <- assoc_results
  d$neglog10 <- -log10(pmax(d$p_adj, .Machine$double.eps))
  ggplot2::ggplot(d, ggplot2::aes(x = factor(.data$component), y = .data$variable,
                                  fill = .data$neglog10)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "#2166ac") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Component", y = "Variable", fill = "-log10(adj p)",
                  title = "Component-metadata associations")
}

.plot_jackstraw_volcano <- function(jackstraw_result, block, component) {
  sig <- rank_features(jackstraw_result = jackstraw_result,
                       mode = "significant",
                       block = block,
                       component = component)
  if (is.null(sig) || nrow(sig) == 0L)
    stop("No jackstraw significance rows available.", call. = FALSE)
  sig$neglog10p <- -log10(pmax(sig$p_value, .Machine$double.eps))
  ggplot2::ggplot(sig, ggplot2::aes(x = .data$feature_index, y = .data$neglog10p,
                                    color = .data$significant)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Feature index", y = "-log10(p)", title = "Jackstraw volcano-style plot")
}

.plot_jackstraw_summary <- function(jackstraw_result) {
  sm <- summarize_components(jackstraw_result = jackstraw_result,
                             summary_type = "significance")
  if (nrow(sm) == 0L) stop("No jackstraw summary available.", call. = FALSE)
  ggplot2::ggplot(sm, ggplot2::aes(x = factor(.data$component), y = .data$n_significant,
                                   fill = factor(.data$block))) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Component", y = "Significant features",
                  fill = "Block", title = "Jackstraw significance summary")
}

.plot_stability_summary <- function(stability_result, style = NULL) {
  sm <- summarize_components(stability_result = stability_result,
                             summary_type = "stability")
  st <- if (is.null(style)) "heatmap" else style

  if (all(c("joint_rank", "frequency") %in% names(sm))) {
    return(
      ggplot2::ggplot(sm, ggplot2::aes(x = factor(.data$joint_rank), y = .data$frequency)) +
        ggplot2::geom_col(fill = "#2166ac") + ggplot2::theme_bw() +
        ggplot2::labs(x = "Joint rank", y = "Frequency", title = "Joint-rank stability")
    )
  }

  if (st == "line") {
    return(
      ggplot2::ggplot(sm, ggplot2::aes(x = .data$component, y = .data$cos_similarity,
                                       color = .data$block, group = .data$block)) +
        ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::theme_bw() +
        ggplot2::labs(x = "Component", y = "Cosine similarity", title = "Loading stability")
    )
  }

  ggplot2::ggplot(sm, ggplot2::aes(x = factor(.data$component), y = .data$block,
                                   fill = .data$cos_similarity)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "#2166ac") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Component", y = "Block", fill = "Cos sim",
                  title = "Stability heatmap")
}


# ---------------------------------------------------------------------------
# associate_components
# ---------------------------------------------------------------------------

#' Associate AJIVE component scores with sample metadata
#'
#' Tests for associations between estimated joint (or individual) component
#' scores from a \code{\link{Rajive}} decomposition and sample-level metadata
#' variables.
#'
#' @section Statistical warnings:
#'
#' \strong{Score estimation error (Finding 4):}
#' Component scores are estimated quantities, not fixed design variables.
#' Estimation error attenuates effect sizes and inflates Type-I error in naive
#' tests. P-values returned by this function do \strong{not} propagate score
#' estimation uncertainty.  Treat results as post-decomposition exploratory
#' associations, not exact fixed-design inference.
#'
#' \strong{Survival split bias (Finding 5):}
#' When \code{mode = "survival"} and a score-split (\code{split = "median"}
#' or \code{split = "tertile"}) is used, the resulting log-rank or Cox p-value
#' is data-adaptive and may be anti-conservative.  Prefer
#' \code{split = "none"} (continuous-score Cox model) for primary inference.
#' Split-based results should be treated as descriptive summaries only.
#'
#' @param ajive_output An object of class \code{"rajive"} (output of
#'   \code{\link{Rajive}}).
#' @param metadata A \code{data.frame} of sample-level metadata; rows must
#'   correspond in order to the samples in \code{ajive_output}.
#' @param scores Optional precomputed scores matrix (samples x components).
#'   When \code{NULL} (default), joint scores are extracted from
#'   \code{ajive_output}.
#' @param mode Character scalar.  One of \code{"continuous"},
#'   \code{"categorical"}, \code{"survival"}, or \code{"batch"}.
#' @param variable Character scalar.  Name of a single column in
#'   \code{metadata} to test.
#' @param variables Character vector.  Names of multiple columns to test.
#'   Ignored when \code{variable} is supplied.
#' @param method Character scalar.  Test method; defaults are Pearson
#'   correlation (\code{"pearson"}) for continuous, Kruskal-Wallis
#'   (\code{"kruskal"}) for categorical, and log-rank (\code{"logrank"})
#'   for survival.
#' @param adjust Character scalar.  P-value adjustment method passed to
#'   \code{\link[stats]{p.adjust}}.  Default \code{"BH"}.
#' @param time_col Character scalar.  Name of the time column in
#'   \code{metadata} (required for \code{mode = "survival"}).
#' @param status_col Character scalar.  Name of the event-status column in
#'   \code{metadata} (required for \code{mode = "survival"}).
#' @param split Character scalar.  Score-split strategy for survival analysis.
#'   One of \code{"none"} (continuous Cox, recommended for primary inference),
#'   \code{"median"}, or \code{"tertile"}.
#' @param type Character scalar.  Source of scores: \code{"joint"} (default)
#'   or \code{"individual"}.
#' @param block Integer or \code{NULL}.  Block index for individual scores.
#' @param component Integer or \code{NULL}.  Restrict analysis to a single
#'   component index.
#' @param ... Reserved for future arguments.
#'
#' @return A \code{data.frame} with one row per variable x component
#'   combination and columns \code{variable}, \code{component}, \code{stat},
#'   \code{p_value}, \code{p_adj}, and \code{method}.
#'
#' @examples
#' \donttest{
#' n   <- 40; pks <- c(30, 20)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#' meta <- data.frame(group = sample(c("A", "B"), n, replace = TRUE))
#' res  <- associate_components(fit, meta, variable = "group",
#'                              mode = "categorical")
#' }
#'
#' @export
associate_components <- function(ajive_output,
                                 metadata,
                                 scores      = NULL,
                                 mode        = c("continuous", "categorical",
                                                  "survival", "batch"),
                                 variable    = NULL,
                                 variables   = NULL,
                                 method      = NULL,
                                 adjust      = "BH",
                                 time_col    = NULL,
                                 status_col  = NULL,
                                 split       = c("none", "median", "tertile"),
                                 type        = c("joint", "individual"),
                                 block       = NULL,
                                 component   = NULL,
                                 ...) {

  mode  <- match.arg(mode)
  type  <- match.arg(type)
  split <- match.arg(split)

  # --- input validation ---
  if (!inherits(ajive_output, "rajive"))
    stop("`ajive_output` must be an object of class \"rajive\".",
         call. = FALSE)
  if (!is.data.frame(metadata))
    stop("`metadata` must be a data.frame.", call. = FALSE)

  # resolve variable list
  vars <- if (!is.null(variable)) variable else variables
  if (is.null(vars) || length(vars) == 0L)
    stop("Supply at least one variable name via `variable` or `variables`.",
         call. = FALSE)
  missing_vars <- setdiff(vars, names(metadata))
  if (length(missing_vars) > 0L)
    stop("Variables not found in `metadata`: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)

  # --- mandatory inferential warning (#4) ---
  message(
    "[associate_components] NOTE: Component scores are estimated quantities. ",
    "Score estimation error is NOT propagated into the returned p-values. ",
    "Treat results as post-decomposition exploratory associations, not exact ",
    "fixed-design inference (StatisticalAudits.md, Finding 4)."
  )

  # --- survival split bias warning (#5) ---
  if (mode == "survival" && split != "none") {
    message(
      "[associate_components] WARNING: Split-based survival inference ",
      "(split = \"", split, "\") is data-adaptive and may be anti-conservative. ",
      "Use split = \"none\" (continuous-score Cox model) for primary inference. ",
      "Split-based results should be labeled as descriptive summaries ",
      "(StatisticalAudits.md, Finding 5)."
    )
  }

  # --- extract scores ---
  if (is.null(scores)) {
    if (type == "joint") {
      sc <- ajive_output$joint_scores
      if (is.null(sc))
        stop("Joint scores not found in `ajive_output`. ",
             "Was `Rajive()` run with a positive joint rank?", call. = FALSE)
    } else {
      if (is.null(block))
        stop("`block` must be specified when `type = \"individual\"`.",
             call. = FALSE)
      sc <- ajive_output$block_decomps[[3L * (as.integer(block) - 1L) + 1L]]$u
      if (is.null(sc))
        stop("Individual scores not found for block ", block, ".", call. = FALSE)
    }
    scores <- sc
  }

  if (is.vector(scores))
    scores <- matrix(scores, ncol = 1L)

  n_comp <- ncol(scores)
  comps  <- if (!is.null(component)) component else seq_len(n_comp)

  if (nrow(scores) != nrow(metadata))
    stop("`metadata` must have the same number of rows as samples in scores (",
         nrow(scores), ").", call. = FALSE)

  # --- run per variable x component ---
  rows <- vector("list", length(vars) * length(comps))
  idx  <- 1L

  for (v in vars) {
    y <- metadata[[v]]
    for (j in comps) {
      x <- scores[, j]
      res_j <- .associate_one(x, y, mode, method, time_col, status_col,
                              split, metadata)
      rows[[idx]] <- data.frame(
        variable  = v,
        component = j,
        stat      = res_j$stat,
        p_value   = res_j$p_value,
        method    = res_j$method,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = adjust)
  out[, c("variable", "component", "stat", "p_value", "p_adj", "method")]
}


# Internal: .associate_one — single variable x component test dispatch
.associate_one <- function(scores_j, y, mode, method, time_col, status_col,
                           split, metadata) {
  switch(mode,

    continuous = {
      meth <- if (!is.null(method)) method else "pearson"
      ct   <- suppressWarnings(
        stats::cor.test(scores_j, as.numeric(y), method = meth)
      )
      list(stat = unname(ct$estimate), p_value = ct$p.value, method = meth)
    },

    categorical = {
      meth <- if (!is.null(method)) method else "kruskal"
      kt   <- stats::kruskal.test(scores_j ~ factor(y))
      list(stat = unname(kt$statistic), p_value = kt$p.value, method = meth)
    },

    survival = {
      if (is.null(time_col) || is.null(status_col))
        stop("`time_col` and `status_col` must be supplied for ",
             "mode = \"survival\".", call. = FALSE)
      # Continuous Cox (primary; split applied below if requested)
      t_vec <- metadata[[time_col]]
      s_vec <- metadata[[status_col]]
      if (!requireNamespace("survival", quietly = TRUE))
        stop("Package 'survival' is required for mode = \"survival\".",
             call. = FALSE)
      s_norm <- if (is.logical(s_vec)) {
        s_vec
      } else if (is.numeric(s_vec)) {
        if (!all(s_vec %in% c(0, 1, NA))) {
          cli::cli_abort(
            c("`status_col` numeric input must contain only 0/1 (or NA).",
              "x" = "Found values outside {0,1}."),
            class = "rajiveplus_invalid_input"
          )
        }
        as.logical(s_vec)
      } else if (is.factor(s_vec) && nlevels(s_vec) == 2L) {
        as.integer(s_vec) == 2L
      } else {
        cli::cli_abort(
          c("`status_col` must be 0/1 numeric, logical, or a 2-level factor.",
            "x" = "Got class {.cls {class(s_vec)}}."),
          class = "rajiveplus_invalid_input"
        )
      }
      surv_obj <- survival::Surv(as.numeric(t_vec), s_norm)
      if (split == "none") {
        fit <- survival::coxph(surv_obj ~ scores_j)
        sm  <- summary(fit)
        list(stat = sm$coefficients[1L, "z"],
             p_value = sm$coefficients[1L, "Pr(>|z|)"],
             method = "cox_continuous")
      } else {
        grp <- if (split == "median") {
          ifelse(scores_j >= stats::median(scores_j), "high", "low")
        } else {
          qs  <- stats::quantile(scores_j, c(1/3, 2/3))
          cut(scores_j, breaks = c(-Inf, qs, Inf),
              labels = c("low", "mid", "high"))
        }
        lr <- survival::survdiff(surv_obj ~ grp)
        p  <- 1 - stats::pchisq(lr$chisq, df = length(unique(grp)) - 1L)
        list(stat = lr$chisq, p_value = p,
             method = paste0("logrank_", split))
      }
    },

    batch = {
      meth <- if (!is.null(method)) method else "kruskal"
      kt   <- stats::kruskal.test(scores_j ~ factor(y))
      list(stat = unname(kt$statistic), p_value = kt$p.value, method = meth)
    }
  )
}


# ---------------------------------------------------------------------------
# assess_stability
# ---------------------------------------------------------------------------

#' Assess stability of AJIVE decomposition via bootstrap resampling
#'
#' Evaluates the stability of the estimated joint rank or block loadings via
#' bootstrap resampling of the sample dimension.
#'
#' @section Procrustes alignment for loading stability (Finding 3):
#' When \code{target = "loadings"}, bootstrap loading matrices are aligned to
#' the reference (full-data) loading matrix using orthogonal Procrustes
#' rotation before any variability summary is computed.  Raw bootstrap loading
#' variation is dominated by rotational indeterminacy and is not interpretable
#' without alignment.  All returned stability summaries are computed
#' \strong{after} Procrustes alignment (StatisticalAudits.md, Finding 3).
#'
#' @param ajive_output An object of class \code{"rajive"} (output of
#'   \code{\link{Rajive}}).  Required for \code{target = "loadings"} and
#'   \code{target = "components"}.
#' @param blocks List of data matrices (same list passed to
#'   \code{\link{Rajive}}).
#' @param initial_signal_ranks Integer vector of initial signal ranks for each
#'   block.
#' @param target Character scalar.  One of \code{"joint_rank"},
#'   \code{"loadings"}, or \code{"components"}.
#' @param method Character scalar.  \code{"bootstrap"} (default) or
#'   \code{"permutation"}.
#' @param B Positive integer.  Number of bootstrap replicates.  Default 100.
#' @param n_perm Positive integer.  Number of permutations (only used when
#'   \code{method = "permutation"}).  Default 100.
#' @param sample_frac Numeric in (0, 1].  Fraction of samples to draw in each
#'   bootstrap replicate.  Default 0.8.
#' @param num_cores Positive integer.  Number of cores for parallel execution.
#'   Default 1 (no parallelism).
#' @param ... Reserved for future arguments.
#'
#' @return A named list whose structure depends on \code{target}:
#' \describe{
#'   \item{\code{"joint_rank"}}{A list with fields \code{rank_distribution}
#'     (integer vector of length \code{B}), \code{rank_table} (frequency
#'     table), and \code{observed_rank} (integer).}
#'   \item{\code{"loadings"}}{A list with one element per block.  Each element
#'     contains \code{mean_loading} (aligned mean loading matrix),
#'     \code{sd_loading} (element-wise standard deviation after alignment), and
#'     \code{cos_similarity} (vector of mean cosine similarities per
#'     component).  All summaries are computed after orthogonal Procrustes
#'     alignment.}
#' }
#'
#' @examples
#' \donttest{
#' n   <- 40; pks <- c(30, 20)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#' stab <- assess_stability(fit, Y$sim_data, c(4, 3),
#'                          target = "joint_rank", B = 50)
#' stab$rank_table
#' }
#'
#' @export
assess_stability <- function(ajive_output = NULL,
                             blocks,
                             initial_signal_ranks,
                             target      = c("joint_rank", "loadings",
                                              "components"),
                             method      = c("bootstrap", "permutation"),
                             B           = 100L,
                             n_perm      = 100L,
                             sample_frac = 0.8,
                             num_cores   = 1L,
                             ...) {

  target <- match.arg(target)
  method <- match.arg(method)

  if (!is.list(blocks) || length(blocks) == 0L)
    stop("`blocks` must be a non-empty list of matrices.", call. = FALSE)
  if (length(initial_signal_ranks) != length(blocks))
    stop("`initial_signal_ranks` must have the same length as `blocks`.",
         call. = FALSE)
  B <- as.integer(B)
  if (B < 2L)
    stop("`B` must be at least 2.", call. = FALSE)
  if (sample_frac <= 0 || sample_frac > 1)
    stop("`sample_frac` must be in (0, 1].", call. = FALSE)

  n_samples <- nrow(blocks[[1L]])

  switch(target,

    joint_rank = {
      rank_draws <- integer(B)
      for (b in seq_len(B)) {
        idx    <- sample.int(n_samples, size = max(2L, floor(sample_frac * n_samples)),
                             replace = TRUE)
        b_list <- lapply(blocks, function(X) X[idx, , drop = FALSE])
        fit_b  <- tryCatch(
          Rajive(b_list, initial_signal_ranks),
          error = function(e) NULL
        )
        rank_draws[b] <- if (is.null(fit_b)) NA_integer_ else fit_b$joint_rank
      }
      list(
        rank_distribution = rank_draws,
        rank_table        = table(rank_draws, useNA = "ifany"),
        observed_rank     = if (!is.null(ajive_output)) ajive_output$joint_rank else NA_integer_
      )
    },

    loadings = {
      if (is.null(ajive_output))
        stop("`ajive_output` must be supplied for target = \"loadings\".",
             call. = FALSE)
      if (!inherits(ajive_output, "rajive"))
        stop("`ajive_output` must be of class \"rajive\".", call. = FALSE)

      K          <- length(blocks)
      joint_rank <- ajive_output$joint_rank
      if (joint_rank == 0L)
        stop("Cannot assess loading stability when joint_rank = 0.", call. = FALSE)

      # Reference loadings: list of K matrices (features x joint_rank)
      # block_decomps layout: individual at 3k-2, joint at 3k-1, noise at 3k
      # $v = features x rank (loadings); $u = samples x rank (scores)
      ref_loadings <- lapply(seq_len(K), function(k) {
        ajive_output$block_decomps[[3L * (k - 1L) + 2L]]$v
      })

      # Bootstrap loading draws: list-of-K, each element is an array
      # (features x joint_rank x B)
      boot_arrays <- lapply(seq_len(K), function(k) {
        array(NA_real_,
              dim = c(nrow(ref_loadings[[k]]), joint_rank, B))
      })

      for (b in seq_len(B)) {
        idx    <- sample.int(n_samples, size = max(2L, floor(sample_frac * n_samples)),
                             replace = TRUE)
        b_list <- lapply(blocks, function(X) X[idx, , drop = FALSE])
        fit_b  <- tryCatch(
          Rajive(b_list, initial_signal_ranks),
          error = function(e) NULL
        )
        if (is.null(fit_b) || fit_b$joint_rank < joint_rank) next

        for (k in seq_len(K)) {
          L_b <- fit_b$block_decomps[[3L * (k - 1L) + 2L]]$v
          # Procrustes-align L_b to reference (Finding #3)
          Q         <- .procrustes_align(ref_loadings[[k]], L_b)
          L_aligned <- L_b %*% Q
          boot_arrays[[k]][, , b] <- L_aligned
        }
      }

      # Summarise across bootstrap replicates
      result <- vector("list", K)
      names(result) <- names(blocks)
      for (k in seq_len(K)) {
        arr         <- boot_arrays[[k]]
        mean_load   <- apply(arr, c(1L, 2L), mean, na.rm = TRUE)
        sd_load     <- apply(arr, c(1L, 2L), stats::sd, na.rm = TRUE)
        ref         <- ref_loadings[[k]]
        cos_sim     <- vapply(seq_len(joint_rank), function(j) {
          ref_j  <- ref[, j]
          mn_j   <- mean_load[, j]
          sum(ref_j * mn_j) /
            (sqrt(sum(ref_j^2)) * sqrt(sum(mn_j^2)) + .Machine$double.eps)
        }, numeric(1L))
        result[[k]] <- list(
          mean_loading  = mean_load,
          sd_loading    = sd_load,
          cos_similarity = cos_sim
        )
      }
      result
    },

    components = {
      if (is.null(ajive_output))
        stop("`ajive_output` must be supplied for target = 'components'.", call. = FALSE)
      if (!inherits(ajive_output, "rajive"))
        stop("`ajive_output` must be of class \"rajive\".", call. = FALSE)

      ref <- ajive_output$joint_scores
      if (is.null(ref) || ncol(ref) == 0L)
        stop("Joint scores are unavailable for component-stability assessment.", call. = FALSE)

      n_comp <- ncol(ref)
      cor_mat <- matrix(NA_real_, nrow = B, ncol = n_comp)

      for (b in seq_len(B)) {
        idx <- sample.int(n_samples, size = max(2L, floor(sample_frac * n_samples)),
                          replace = TRUE)
        b_list <- lapply(blocks, function(X) X[idx, , drop = FALSE])
        fit_b <- tryCatch(
          Rajive(b_list, initial_signal_ranks),
          error = function(e) NULL
        )
        if (is.null(fit_b) || is.null(fit_b$joint_scores)) next
        bs <- fit_b$joint_scores
        n_use <- min(ncol(bs), n_comp)
        if (n_use < 1L) next
        # Align bootstrap components to the reference subspace to resolve
        # sign and order ambiguity before per-component correlations.
        ref_sub <- ref[idx, seq_len(n_use), drop = FALSE]
        bs_sub  <- bs[, seq_len(n_use), drop = FALSE]
        Q <- .procrustes_align(ref_sub, bs_sub)
        bs_aligned <- bs_sub %*% Q
        for (j in seq_len(n_use)) {
          # Keep absolute correlation for residual sign indeterminacy.
          cor_mat[b, j] <- suppressWarnings(
            abs(stats::cor(ref_sub[, j], bs_aligned[, j],
                           use = "pairwise.complete.obs"))
          )
        }
      }

      data.frame(
        component = seq_len(n_comp),
        mean_correlation = apply(cor_mat, 2L, mean, na.rm = TRUE),
        sd_correlation = apply(cor_mat, 2L, stats::sd, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )
}


# ---------------------------------------------------------------------------
# Internal: orthogonal Procrustes alignment
# ---------------------------------------------------------------------------

# .procrustes_align(A, B) returns an orthogonal rotation matrix Q such that
# B %*% Q is as close as possible to A in Frobenius norm.
# Solution: Q = V %*% t(U) where svd(t(A) %*% B) = U D V'.
.procrustes_align <- function(A, B) {
  P     <- crossprod(A, B)          # t(A) %*% B
  sv    <- svd(P)
  sv$v %*% t(sv$u)
}


# ---------------------------------------------------------------------------
# rank_features
# ---------------------------------------------------------------------------

#' Rank and summarize features from a Rajive decomposition
#'
#' Returns a tidy \code{data.frame} of features ranked by their loading,
#' proportional contribution, pairwise feature-set overlap, or jackstraw
#' significance, depending on \code{mode}.
#'
#' @param ajive_output An object of class \code{"rajive"} (output of
#'   \code{\link{Rajive}}).  Required for \code{mode} in
#'   \code{"top_loadings"}, \code{"contribution"}, and \code{"overlap"}.
#'   May be \code{NULL} for \code{mode = "significant"}.
#' @param blocks Not currently used; reserved for future feature-name lookup
#'   against the original data matrices.
#' @param jackstraw_result An object of class \code{"jackstraw_rajive"} (output
#'   of \code{\link{jackstraw_rajive}}).  Required for
#'   \code{mode = "significant"}.
#' @param mode Character scalar.  One of:
#'   \describe{
#'     \item{\code{"top_loadings"}}{Return the top \code{n} features per
#'       block \eqn{\times} component by absolute loading.}
#'     \item{\code{"contribution"}}{Return the top \code{n} features per
#'       block \eqn{\times} component by their normalized contribution score.
#'       For \code{method = "abs_loading"}: \eqn{|w_i| / \sum |w_j|}.
#'       For \code{method = "variance_contrib"}:
#'       \eqn{w_i^2 / \sum w_j^2}.}
#'     \item{\code{"overlap"}}{Compute pairwise Jaccard or overlap-coefficient
#'       between the top-\code{n} feature sets of each pair of blocks for each
#'       requested component.}
#'     \item{\code{"significant"}}{Return all features from a jackstraw result
#'       with their p-values and significance calls, sorted by block, component,
#'       and p-value.}
#'   }
#' @param type Character.  \code{"joint"} (default) or \code{"individual"}.
#'   Selects which loading matrix is used.
#' @param block Integer or \code{NULL}.  Block index to restrict results to.
#'   When \code{NULL} (default), all blocks are included.
#' @param component Integer or \code{NULL}.  Component index to restrict
#'   results to.  When \code{NULL} (default), all components are included.
#' @param n Integer.  Number of top features returned per block
#'   \eqn{\times} component (for \code{mode} in \code{"top_loadings"},
#'   \code{"contribution"}).  Also the feature-set size used in
#'   \code{mode = "overlap"}.  Default 20.
#' @param signed Logical.  When \code{TRUE} (default), the \code{loading}
#'   column retains the sign of the original loading.  When \code{FALSE},
#'   the absolute value is stored.  Ignored for \code{mode} in
#'   \code{"overlap"} and \code{"significant"}.
#' @param method Character or \code{NULL}.  Scoring or overlap method.
#'   For \code{mode} in \code{"top_loadings"} and \code{"contribution"}:
#'   \code{"abs_loading"} (default) or \code{"variance_contrib"}.
#'   For \code{mode = "overlap"}: \code{"jaccard"} (default) or
#'   \code{"overlap_coef"}.
#'   Ignored for \code{mode = "significant"}.
#' @param top_n Optional positive integer.  When supplied to
#'   \code{summarize_significant_vars}, retains only the top \code{top_n}
#'   features (by smallest p-value) per block \eqn{\times} component group.
#'   Ignored by \code{rank_features} itself and other aliases.
#' @param ... Reserved for future arguments.
#'
#' @return A \code{data.frame}.  Column set depends on \code{mode}:
#'   \describe{
#'     \item{\code{"top_loadings"}}{
#'       \code{block}, \code{component}, \code{type},
#'       \code{feature_index}, \code{feature_name},
#'       \code{loading}, \code{abs_loading}, \code{rank}.}
#'     \item{\code{"contribution"}}{
#'       Same as \code{"top_loadings"} plus \code{contribution}.}
#'     \item{\code{"overlap"}}{
#'       \code{component}, \code{type}, \code{block_i}, \code{block_j},
#'       \code{top_n}, \code{method}, \code{n_intersect},
#'       \code{overlap_score}.}
#'     \item{\code{"significant"}}{
#'       \code{block}, \code{component}, \code{feature_index},
#'       \code{p_value}, \code{p_adj}, \code{significant}.}
#'   }
#'   Returns \code{NULL} (with a message) when no valid block/component
#'   combinations are found.
#'
#' @examples
#' \donttest{
#' n <- 40; pks <- c(30, 20)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#'
#' # Top 10 joint loadings across all blocks
#' rank_features(fit, mode = "top_loadings", type = "joint", n = 10)
#'
#' # Proportional contribution (variance) for individual components
#' rank_features(fit, mode = "contribution", type = "individual",
#'               block = 1, component = 1, method = "variance_contrib")
#'
#' # Pairwise feature-set overlap across all blocks for joint component 1
#' rank_features(fit, mode = "overlap", type = "joint", component = 1, n = 15)
#'
#' # Significant jackstraw features
#' jsr <- jackstraw_rajive(fit, Y$sim_data)
#' rank_features(jackstraw_result = jsr, mode = "significant", block = 1)
#' }
#'
#' @export
rank_features <- function(ajive_output    = NULL,
                          blocks          = NULL,
                          jackstraw_result = NULL,
                          mode      = c("top_loadings", "contribution",
                                        "overlap", "significant"),
                          type      = c("joint", "individual"),
                          block     = NULL,
                          component = NULL,
                          n         = 20L,
                          signed    = TRUE,
                          method    = NULL,
                          ...) {
  mode <- match.arg(mode)
  type <- match.arg(type)
  n    <- as.integer(n)
  if (n < 1L) stop("`n` must be a positive integer.", call. = FALSE)

  method <- .rf_resolve_method(method, mode)

  switch(mode,
    top_loadings = {
      if (!inherits(ajive_output, "rajive"))
        stop("`ajive_output` must be an object of class \"rajive\".", call. = FALSE)
      .rf_loading_df(ajive_output, type, block, component, n, signed, method,
                     add_contribution = FALSE)
    },
    contribution = {
      if (!inherits(ajive_output, "rajive"))
        stop("`ajive_output` must be an object of class \"rajive\".", call. = FALSE)
      .rf_loading_df(ajive_output, type, block, component, n, signed, method,
                     add_contribution = TRUE)
    },
    overlap = {
      if (!inherits(ajive_output, "rajive"))
        stop("`ajive_output` must be an object of class \"rajive\".", call. = FALSE)
      .rf_overlap(ajive_output, type, block, component, n, method)
    },
    significant = {
      .rf_significant(jackstraw_result, block, component)
    }
  )
}


# ---------------------------------------------------------------------------
# summarize_components
# ---------------------------------------------------------------------------

#' Summarize component-level outputs from rajiveplus helper workflows
#'
#' Aggregates common outputs into compact tables suitable for reporting and
#' downstream plotting. Supported summaries include decomposition structure,
#' variance-explained tables, jackstraw significance counts, metadata
#' associations, stability outputs, and AJIVE rank diagnostics.
#'
#' @param ajive_output Optional object of class \code{"rajive"}.
#' @param jackstraw_result Optional object of class \code{"jackstraw_rajive"}.
#' @param assoc_results Optional association result table from
#'   \code{\link{associate_components}}.
#' @param stability_result Optional stability result object from
#'   \code{\link{assess_stability}}.
#' @param blocks Optional list of block matrices (needed for
#'   \code{summary_type = "variance"}).
#' @param summary_type Character scalar. One of \code{"components"},
#'   \code{"variance"}, \code{"significance"}, \code{"associations"},
#'   \code{"stability"}, or \code{"diagnostics"}.
#' @param block Optional integer block index filter.
#' @param component Optional integer component index filter.
#' @param top_n Optional integer to keep only top rows per summary.
#' @param ... Reserved for future extensions.
#'
#' @return A \code{data.frame} with stable columns for the requested summary
#'   type.
#'
#' @examples
#' \donttest{
#' n <- 40; pks <- c(30, 20)
#' Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                     n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#' summarize_components(ajive_output = fit, summary_type = "components")
#' summarize_components(ajive_output = fit, summary_type = "diagnostics")
#' }
#'
#' @export
summarize_components <- function(ajive_output = NULL,
                                 jackstraw_result = NULL,
                                 assoc_results = NULL,
                                 stability_result = NULL,
                                 blocks = NULL,
                                 summary_type = c("components", "variance",
                                                  "significance", "associations",
                                                  "stability", "diagnostics"),
                                 block = NULL,
                                 component = NULL,
                                 top_n = NULL,
                                 ...) {
  summary_type <- match.arg(summary_type)

  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    if (is.na(top_n) || top_n < 1L)
      stop("`top_n` must be a positive integer when supplied.", call. = FALSE)
  }

  out <- switch(summary_type,

    components = {
      if (!inherits(ajive_output, "rajive"))
        stop("`ajive_output` must be an object of class \"rajive\" for summary_type = 'components'.",
             call. = FALSE)

      K <- length(ajive_output$block_decomps) %/% 3L
      blocks_k <- if (is.null(block)) seq_len(K) else as.integer(block)
      rows <- list()

      for (k in blocks_k) {
        if (k < 1L || k > K) next
        ind_v <- ajive_output$block_decomps[[3L * (k - 1L) + 1L]]$v
        jnt_v <- ajive_output$block_decomps[[3L * (k - 1L) + 2L]]$v
        ind_r <- if (is.null(ind_v)) 0L else ncol(ind_v)
        jnt_r <- if (is.null(jnt_v)) 0L else ncol(jnt_v)

        add_rows <- function(tp, rank_k, vmat) {
          if (rank_k < 1L) return(NULL)
          comps <- if (is.null(component)) seq_len(rank_k) else as.integer(component)
          comps <- comps[comps >= 1L & comps <= rank_k]
          if (length(comps) == 0L) return(NULL)
          data.frame(
            block = k,
            type = tp,
            component = comps,
            block_rank = rank_k,
            joint_rank = as.integer(ajive_output$joint_rank),
            n_features = nrow(vmat),
            stringsAsFactors = FALSE
          )
        }

        rows[[length(rows) + 1L]] <- add_rows("joint", jnt_r, jnt_v)
        rows[[length(rows) + 1L]] <- add_rows("individual", ind_r, ind_v)
      }
      rows <- Filter(Negate(is.null), rows)
      if (length(rows) == 0L) {
        data.frame(block = integer(0), type = character(0), component = integer(0),
                   block_rank = integer(0), joint_rank = integer(0),
                   n_features = integer(0), stringsAsFactors = FALSE)
      } else {
        do.call(rbind, rows)
      }
    },

    variance = {
      if (!inherits(ajive_output, "rajive"))
        stop("`ajive_output` must be an object of class \"rajive\" for summary_type = 'variance'.",
             call. = FALSE)
      if (is.null(blocks) || !is.list(blocks))
        stop("`blocks` must be supplied as a list for summary_type = 'variance'.",
             call. = FALSE)

      vt <- as.data.frame(showVarExplained_robust(ajive_output, blocks),
                          stringsAsFactors = FALSE)
      names(vt) <- make.names(names(vt), unique = TRUE)
      if (!is.null(block)) {
        b <- as.integer(block)
        if ("Block" %in% names(vt)) {
          vt <- vt[vt$Block %in% b | vt$Block %in% paste0("block", b), , drop = FALSE]
        }
      }
      vt
    },

    significance = {
      sig <- rank_features(jackstraw_result = jackstraw_result,
                           mode = "significant",
                           block = block,
                           component = component)
      if (is.null(sig) || nrow(sig) == 0L) {
        data.frame(block = integer(0), component = integer(0), n_features = integer(0),
                   n_significant = integer(0), prop_significant = numeric(0),
                   min_p_value = numeric(0), min_p_adj = numeric(0),
                   stringsAsFactors = FALSE)
      } else {
        grp <- split(sig, paste(sig$block, sig$component))
        out_s <- do.call(rbind, lapply(grp, function(df) {
          data.frame(
            block = df$block[1L],
            component = df$component[1L],
            n_features = nrow(df),
            n_significant = sum(df$significant, na.rm = TRUE),
            prop_significant = mean(df$significant, na.rm = TRUE),
            min_p_value = min(df$p_value, na.rm = TRUE),
            min_p_adj = min(df$p_adj, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }))
        out_s[order(out_s$block, out_s$component), , drop = FALSE]
      }
    },

    associations = {
      if (is.null(assoc_results) || !is.data.frame(assoc_results))
        stop("`assoc_results` must be a data.frame for summary_type = 'associations'.",
             call. = FALSE)
      req <- c("variable", "component", "p_value", "p_adj")
      missing_req <- setdiff(req, names(assoc_results))
      if (length(missing_req) > 0L)
        stop("`assoc_results` is missing required columns: ",
             paste(missing_req, collapse = ", "), call. = FALSE)

      a <- assoc_results
      if (!is.null(component)) a <- a[a$component %in% as.integer(component), , drop = FALSE]
      if (!is.null(block) && "block" %in% names(a)) {
        a <- a[a$block %in% as.integer(block), , drop = FALSE]
      }
      if (nrow(a) == 0L) {
        data.frame(variable = character(0), component = integer(0), n_tests = integer(0),
                   min_p_value = numeric(0), min_p_adj = numeric(0), method = character(0),
                   stringsAsFactors = FALSE)
      } else {
        grp <- split(a, paste(a$variable, a$component, sep = "||"))
        do.call(rbind, lapply(grp, function(df) {
          data.frame(
            variable = df$variable[1L],
            component = df$component[1L],
            n_tests = nrow(df),
            min_p_value = min(df$p_value, na.rm = TRUE),
            min_p_adj = min(df$p_adj, na.rm = TRUE),
            method = paste(unique(df$method), collapse = ";"),
            stringsAsFactors = FALSE
          )
        }))
      }
    },

    stability = {
      if (is.null(stability_result) || !is.list(stability_result))
        stop("`stability_result` must be a list from assess_stability() for summary_type = 'stability'.",
             call. = FALSE)

      if (all(c("rank_distribution", "rank_table") %in% names(stability_result))) {
        rt <- as.data.frame(stability_result$rank_table, stringsAsFactors = FALSE)
        if (ncol(rt) < 2L) {
          data.frame(joint_rank = integer(0), frequency = integer(0), proportion = numeric(0),
                     stringsAsFactors = FALSE)
        } else {
          names(rt)[1:2] <- c("joint_rank", "frequency")
          rt$joint_rank <- suppressWarnings(as.integer(as.character(rt$joint_rank)))
          rt$proportion <- rt$frequency / sum(rt$frequency)
          rt[order(rt$joint_rank), , drop = FALSE]
        }
      } else {
        block_names <- names(stability_result)
        if (is.null(block_names) || any(block_names == ""))
          block_names <- paste0("block", seq_along(stability_result))
        rows <- list()
        for (k in seq_along(stability_result)) {
          sk <- stability_result[[k]]
          if (is.null(sk$cos_similarity)) next
          cs <- as.numeric(sk$cos_similarity)
          comps <- seq_along(cs)
          if (!is.null(component)) {
            keep <- comps %in% as.integer(component)
            cs <- cs[keep]
            comps <- comps[keep]
          }
          if (length(cs) == 0L) next
          rows[[length(rows) + 1L]] <- data.frame(
            block = block_names[k],
            component = comps,
            cos_similarity = cs,
            stringsAsFactors = FALSE
          )
        }
        if (length(rows) == 0L) {
          data.frame(block = character(0), component = integer(0),
                     cos_similarity = numeric(0), stringsAsFactors = FALSE)
        } else {
          do.call(rbind, rows)
        }
      }
    },

    diagnostics = {
      if (!inherits(ajive_output, "rajive"))
        stop("`ajive_output` must be an object of class \"rajive\" for summary_type = 'diagnostics'.",
             call. = FALSE)
      if (!is.null(block))
        warning("`block` is ignored for diagnostics summaries.", call. = FALSE)

      d <- extract_components(ajive_output, what = "rank_diagnostics", format = "long")
      if (!is.null(component)) d <- d[d$component_index %in% as.integer(component), , drop = FALSE]
      d
    }
  )

  if (!is.null(top_n) && is.data.frame(out) && nrow(out) > 0L) {
    if (all(c("block", "component") %in% names(out))) {
      split_key <- paste(out$block, out$component, sep = "||")
      out <- do.call(rbind, lapply(split(out, split_key), function(df) {
        df[seq_len(min(nrow(df), top_n)), , drop = FALSE]
      }))
      rownames(out) <- NULL
    } else {
      out <- out[seq_len(min(nrow(out), top_n)), , drop = FALSE]
    }
  }

  out
}


# ---------------------------------------------------------------------------
# export_results
# ---------------------------------------------------------------------------

#' Export summarized rajiveplus outputs to disk
#'
#' Writes a data frame, list, or plot object produced by the rajiveplus
#' analysis pipeline to a file on disk.  Plain-text formats (
#' \code{"csv"}, \code{"tsv"}, \code{"md"}, \code{"html"}) require
#' \code{x} to be a \code{data.frame}.  For arbitrary R objects use
#' \code{format = "rds"}.  When \code{include_plots = TRUE} and \code{x} is a
#' list containing \code{ggplot} objects, each plot is also saved as a PNG file
#' in the same directory as \code{path}.
#'
#' @param x Object to export. Typically a data frame, list, or plot object.
#' @param path Character. Output file path (directory is created if it does not
#'   exist).
#' @param format Export format: \code{"csv"}, \code{"tsv"}, \code{"html"},
#'   \code{"md"}, or \code{"rds"}.
#' @param report Character. Report label embedded in the file header for
#'   bookkeeping purposes.
#' @param include_plots Logical; if \code{TRUE} and \code{x} is a list
#'   containing ggplot objects, write PNG files next to the main output.
#' @param ... Reserved for future extensions.
#'
#' @return Invisibly returns \code{path}.
#'
#' @examples
#' \donttest{
#' n   <- 40; pks <- c(30, 20)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#' feat <- rank_features(fit, mode = "top_loadings", type = "joint", n = 5L)
#' tmp  <- tempfile(fileext = ".csv")
#' export_results(feat, path = tmp, format = "csv")
#' }
#'
#' @export
export_results <- function(x,
                           path,
                           format = c("csv", "tsv", "html", "md", "rds"),
                           report = c("top_features", "component_summary",
                                      "associations", "full"),
                           include_plots = TRUE,
                           ...) {
  format <- match.arg(format)
  report <- match.arg(report)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (format == "rds") {
    saveRDS(x, path)
  } else if (format %in% c("csv", "tsv")) {
    if (!is.data.frame(x))
      stop("For csv/tsv export, `x` must be a data.frame.", call. = FALSE)
    sep <- if (format == "csv") "," else "\t"
    utils::write.table(x, file = path, sep = sep, quote = TRUE,
                       row.names = FALSE, col.names = TRUE)
  } else if (format == "md") {
    if (is.data.frame(x)) {
      txt <- paste(capture.output(print(x)), collapse = "\n")
      writeLines(c(paste0("# ", report), "", "```", txt, "```"), con = path)
    } else {
      writeLines(c(paste0("# ", report), "", "```", capture.output(str(x)), "```"), con = path)
    }
  } else if (format == "html") {
    body <- if (is.data.frame(x)) {
      paste("<pre>", paste(capture.output(print(x)), collapse = "\n"), "</pre>")
    } else {
      paste("<pre>", paste(capture.output(str(x)), collapse = "\n"), "</pre>")
    }
    html <- c("<html><head><meta charset='utf-8'></head><body>",
              paste0("<h1>", report, "</h1>"), body, "</body></html>")
    writeLines(html, con = path)
  }

  if (isTRUE(include_plots) && is.list(x)) {
    gg_idx <- which(vapply(x, inherits, logical(1L), what = "ggplot"))
    if (length(gg_idx) > 0L) {
      base <- sub("\\.[^.]+$", "", basename(path))
      out_dir <- dirname(path)
      for (i in seq_along(gg_idx)) {
        nm <- names(x)[gg_idx[i]]
        if (is.null(nm) || identical(nm, "")) nm <- paste0("plot", i)
        png_file <- file.path(out_dir, paste0(base, "_", nm, ".png"))
        ggplot2::ggsave(filename = png_file, plot = x[[gg_idx[i]]],
                        width = 7, height = 5, dpi = 120)
      }
    }
  }

  invisible(path)
}


# ---------------------------------------------------------------------------
# rajive_report
# ---------------------------------------------------------------------------

#' Build a one-shot interpretation report from Rajive outputs
#'
#' Assembles an HTML or Markdown interpretation report by calling a set of
#' summary and visualization functions on a \code{\link{Rajive}} result and,
#' optionally, a \code{\link{jackstraw_rajive}} significance result and a
#' sample metadata frame.  Each requested \code{section} is written in turn;
#' sections are silently skipped when the required inputs are absent.  The
#' output file is written to \code{output_file}.
#'
#' @param ajive_output An object of class \code{"rajive"} (output of
#'   \code{\link{Rajive}}).
#' @param blocks List of original block matrices passed to variance and
#'   feature summary sections.
#' @param metadata Optional \code{data.frame} of sample-level metadata.  When
#'   supplied, the associations section is populated for the first column.
#' @param jackstraw_result Optional object of class \code{"jackstraw_rajive"}
#'   (output of \code{\link{jackstraw_rajive}}).  When supplied, a
#'   significance summary section is included in the report.
#' @param output_file Character. Output file path.  The extension determines
#'   whether an HTML (\.html) or Markdown (\.md) report is written.  Defaults
#'   to \code{"rajive_interpretation_report.html"}.
#' @param sections Character vector of sections to include in the report.
#'   Supported values: \code{"overview"}, \code{"variance"},
#'   \code{"features"}, \code{"associations"}, \code{"stability"}.
#' @param ... Reserved for future extensions.
#'
#' @return Invisibly returns \code{output_file}.
#'
#' @examples
#' \donttest{
#' n   <- 40; pks <- c(30, 20)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#' tmp <- tempfile(fileext = ".html")
#' rajive_report(fit, blocks = Y$sim_data, output_file = tmp,
#'               sections = c("overview", "variance", "features"))
#' }
#'
#' @export
rajive_report <- function(ajive_output,
                          blocks,
                          metadata = NULL,
                          jackstraw_result = NULL,
                          output_file = "rajive_interpretation_report.html",
                          sections = c("overview", "variance", "features",
                                       "associations", "stability"),
                          ...) {
  if (!inherits(ajive_output, "rajive"))
    stop("`ajive_output` must be an object of class 'rajive'.", call. = FALSE)
  if (!is.list(blocks))
    stop("`blocks` must be a list.", call. = FALSE)

  lines <- c("# rajiveplus Interpretation Report", "")

  if ("overview" %in% sections) {
    comp <- summarize_components(ajive_output = ajive_output,
                                 summary_type = "components")
    lines <- c(lines, "## Overview", "", paste(capture.output(print(comp)), collapse = "\n"), "")
  }

  if ("variance" %in% sections) {
    var_tbl <- summarize_components(ajive_output = ajive_output,
                                    blocks = blocks,
                                    summary_type = "variance")
    lines <- c(lines, "## Variance", "", paste(capture.output(print(var_tbl)), collapse = "\n"), "")
  }

  if ("features" %in% sections) {
    feat <- rank_features(ajive_output = ajive_output, mode = "top_loadings",
                          type = "joint", n = 10L)
    lines <- c(lines, "## Top Features", "", paste(capture.output(print(feat)), collapse = "\n"), "")
  }

  if ("associations" %in% sections && !is.null(metadata) && is.data.frame(metadata)) {
    vars <- names(metadata)
    if (length(vars) > 0L) {
      assoc <- tryCatch(
        associate_components(ajive_output = ajive_output,
                             metadata = metadata,
                             variable = vars[1L],
                             mode = if (is.numeric(metadata[[vars[1L]]])) "continuous" else "categorical"),
        error = function(e) NULL
      )
      if (!is.null(assoc)) {
        lines <- c(lines, "## Associations", "", paste(capture.output(print(assoc)), collapse = "\n"), "")
      }
    }
  }

  if ("stability" %in% sections) {
    lines <- c(lines, "## Stability", "", "Run assess_stability() with desired bootstrap settings for full stability summaries.", "")
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.html$", output_file, ignore.case = TRUE)) {
    html <- c("<html><head><meta charset='utf-8'></head><body>",
              "<pre>", lines, "</pre>", "</body></html>")
    writeLines(html, con = output_file)
  } else {
    writeLines(lines, con = output_file)
  }

  invisible(output_file)
}


# ---------------------------------------------------------------------------
# Compatibility aliases
# ---------------------------------------------------------------------------

#' @rdname associate_components
#' @export
associate_scores_continuous <- function(ajive_output, metadata, variable,
                                        type = c("joint", "individual"),
                                        block = NULL, component = NULL,
                                        adjust = "BH") {
  associate_components(ajive_output = ajive_output,
                       metadata = metadata,
                       variable = variable,
                       mode = "continuous",
                       type = match.arg(type),
                       block = block,
                       component = component,
                       adjust = adjust)
}

#' @rdname associate_components
#' @export
associate_scores_categorical <- function(ajive_output, metadata, variable,
                                         type = c("joint", "individual"),
                                         block = NULL, component = NULL,
                                         adjust = "BH") {
  associate_components(ajive_output = ajive_output,
                       metadata = metadata,
                       variable = variable,
                       mode = "categorical",
                       type = match.arg(type),
                       block = block,
                       component = component,
                       adjust = adjust)
}

#' @rdname associate_components
#' @export
associate_scores_survival <- function(ajive_output, metadata,
                                      time_col, status_col,
                                      split = c("none", "median", "tertile"),
                                      type = c("joint", "individual"),
                                      block = NULL, component = NULL,
                                      adjust = "BH") {
  associate_components(ajive_output = ajive_output,
                       metadata = metadata,
                       variable = time_col,
                       mode = "survival",
                       time_col = time_col,
                       status_col = status_col,
                       split = match.arg(split),
                       type = match.arg(type),
                       block = block,
                       component = component,
                       adjust = adjust)
}

#' @rdname assess_stability
#' @export
bootstrap_joint_rank <- function(ajive_output = NULL,
                                 blocks,
                                 initial_signal_ranks,
                                 B = 100L,
                                 sample_frac = 0.8,
                                 num_cores = 1L) {
  assess_stability(ajive_output = ajive_output,
                   blocks = blocks,
                   initial_signal_ranks = initial_signal_ranks,
                   target = "joint_rank",
                   B = B,
                   sample_frac = sample_frac,
                   num_cores = num_cores)
}

#' @rdname assess_stability
#' @export
bootstrap_loading_stability <- function(ajive_output,
                                        blocks,
                                        initial_signal_ranks,
                                        B = 100L,
                                        sample_frac = 0.8,
                                        num_cores = 1L) {
  assess_stability(ajive_output = ajive_output,
                   blocks = blocks,
                   initial_signal_ranks = initial_signal_ranks,
                   target = "loadings",
                   B = B,
                   sample_frac = sample_frac,
                   num_cores = num_cores)
}

#' @rdname plot_components
#' @export
plot_stability_heatmap <- function(stability_result, style = "heatmap", ...) {
  plot_components(stability_result = stability_result,
                  plot_type = "stability",
                  style = style,
                  ...)
}


# ---------------------------------------------------------------------------
# Aliases (thin wrappers)
# ---------------------------------------------------------------------------

#' @rdname rank_features
#' @export
get_top_loadings <- function(ajive_output, block = NULL, component = NULL,
                             type = c("joint", "individual"),
                             n = 20L, signed = TRUE) {
  rank_features(ajive_output, mode = "top_loadings", type = match.arg(type),
                block = block, component = component, n = n, signed = signed)
}

#' @rdname rank_features
#' @export
get_feature_contributions <- function(ajive_output, block = NULL, component = NULL,
                                      type   = c("joint", "individual"),
                                      method = c("abs_loading",
                                                 "variance_contrib")) {
  rank_features(ajive_output, mode = "contribution", type = match.arg(type),
                block = block, component = component,
                method = match.arg(method))
}

#' @rdname rank_features
#' @export
compare_feature_sets_across_blocks <- function(ajive_output, component = NULL,
                                               type   = c("joint", "individual"),
                                               n      = 50L,
                                               method = c("jaccard",
                                                          "overlap_coef")) {
  rank_features(ajive_output, mode = "overlap", type = match.arg(type),
                component = component, n = n, method = match.arg(method))
}

#' @rdname rank_features
#' @export
summarize_significant_vars <- function(jackstraw_result, block = NULL,
                                       component = NULL, top_n = NULL) {
  out <- rank_features(jackstraw_result = jackstraw_result,
                       mode = "significant", block = block,
                       component = component)
  if (!is.null(top_n) && !is.null(out)) {
    grp <- paste(out$block, out$component)
    out <- do.call(rbind, lapply(split(out, grp), function(df) {
      df <- df[order(df$p_value), ]
      head(df, as.integer(top_n))
    }))
    rownames(out) <- NULL
  }
  out
}


# ---------------------------------------------------------------------------
# Internal: rank_features workers
# ---------------------------------------------------------------------------

.rf_resolve_method <- function(method, mode) {
  if (is.null(method)) {
    return(switch(mode,
      top_loadings = "abs_loading",
      contribution = "abs_loading",
      overlap      = "jaccard",
      significant  = NULL
    ))
  }
  method <- match.arg(method, c("abs_loading", "variance_contrib",
                                "jaccard", "overlap_coef"))
  if (mode %in% c("top_loadings", "contribution") &&
      method %in% c("jaccard", "overlap_coef"))
    stop(sprintf(
      "method = '%s' is not valid for mode = '%s'. Use 'abs_loading' or 'variance_contrib'.",
      method, mode), call. = FALSE)
  if (mode == "overlap" && method %in% c("abs_loading", "variance_contrib"))
    stop(sprintf(
      "method = '%s' is not valid for mode = 'overlap'. Use 'jaccard' or 'overlap_coef'.",
      method), call. = FALSE)
  method
}

# Returns loading matrix (features x components) for block k.
# Individual: index 3*(k-1)+1; joint: index 3*(k-1)+2.
.rf_get_loadings <- function(ajive_output, k, type) {
  k   <- as.integer(k)
  idx <- if (type == "individual") 3L * (k - 1L) + 1L else 3L * (k - 1L) + 2L
  bd  <- ajive_output$block_decomps[[idx]]
  if (is.null(bd) || is.null(bd$v) || ncol(bd$v) == 0L) return(NULL)
  bd$v
}

.rf_loading_df <- function(ajive_output, type, block, component, n,
                           signed, method, add_contribution) {
  K        <- length(ajive_output$block_decomps) %/% 3L
  blocks_k <- if (is.null(block)) seq_len(K) else as.integer(block)

  rows <- list()
  for (k in blocks_k) {
    mat <- .rf_get_loadings(ajive_output, k, type)
    if (is.null(mat)) next
    n_comp     <- ncol(mat)
    feat_names <- rownames(mat)
    comps      <- if (is.null(component)) seq_len(n_comp) else as.integer(component)

    for (j in comps) {
      if (j > n_comp) next
      load_j  <- mat[, j]
      score_j <- if (method == "variance_contrib") {
        ss <- sum(load_j^2)
        if (ss == 0) rep(0, length(load_j)) else load_j^2 / ss
      } else {
        s1 <- sum(abs(load_j))
        if (s1 == 0) rep(0, length(load_j)) else abs(load_j) / s1
      }
      n_take <- min(n, length(load_j))
      ord    <- order(score_j, decreasing = TRUE)[seq_len(n_take)]

      row_df <- data.frame(
        block         = k,
        component     = j,
        type          = type,
        feature_index = ord,
        feature_name  = if (!is.null(feat_names)) feat_names[ord] else NA_character_,
        loading       = if (signed) load_j[ord] else abs(load_j[ord]),
        abs_loading   = abs(load_j[ord]),
        stringsAsFactors = FALSE
      )
      if (add_contribution) row_df$contribution <- score_j[ord]
      row_df$rank <- seq_len(n_take)
      rows[[length(rows) + 1L]] <- row_df
    }
  }

  if (length(rows) == 0L) {
    message("No valid block/component combinations found. Returning NULL.")
    return(NULL)
  }
  do.call(rbind, rows)
}

.rf_overlap <- function(ajive_output, type, block, component, n, method) {
  K        <- length(ajive_output$block_decomps) %/% 3L
  blocks_k <- if (is.null(block)) seq_len(K) else as.integer(block)
  if (length(blocks_k) < 2L)
    stop("mode = 'overlap' requires at least two blocks.", call. = FALSE)

  # use first valid block to infer component count
  ref_mat <- NULL
  for (k in blocks_k) {
    ref_mat <- .rf_get_loadings(ajive_output, k, type)
    if (!is.null(ref_mat)) break
  }
  if (is.null(ref_mat))
    stop("No valid loading matrices found for the requested blocks.", call. = FALSE)
  n_comp <- ncol(ref_mat)
  comps  <- if (is.null(component)) seq_len(n_comp) else as.integer(component)

  rows <- list()
  for (j in comps) {
    feat_sets <- lapply(blocks_k, function(k) {
      mat <- .rf_get_loadings(ajive_output, k, type)
      if (is.null(mat) || j > ncol(mat)) return(integer(0L))
      n_take <- min(n, nrow(mat))
      order(abs(mat[, j]), decreasing = TRUE)[seq_len(n_take)]
    })

    nb <- length(blocks_k)
    for (i in seq_len(nb - 1L)) {
      for (jj in (i + 1L):nb) {
        A       <- feat_sets[[i]]
        B       <- feat_sets[[jj]]
        n_inter <- length(intersect(A, B))
        n_union <- length(union(A, B))
        score   <- if (method == "jaccard") {
          if (n_union == 0L) NA_real_ else n_inter / n_union
        } else {
          denom <- min(length(A), length(B))
          if (denom == 0L) NA_real_ else n_inter / denom
        }
        rows[[length(rows) + 1L]] <- data.frame(
          component     = j,
          type          = type,
          block_i       = blocks_k[i],
          block_j       = blocks_k[jj],
          top_n         = n,
          method        = method,
          n_intersect   = n_inter,
          overlap_score = score,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

.rf_significant <- function(jackstraw_result, block, component) {
  if (is.null(jackstraw_result))
    stop("`jackstraw_result` must be provided for mode = \"significant\".",
         call. = FALSE)
  if (!inherits(jackstraw_result, "jackstraw_rajive"))
    stop("`jackstraw_result` must be of class \"jackstraw_rajive\".",
         call. = FALSE)

  block_names <- names(jackstraw_result)
  use_blocks  <- if (!is.null(block)) {
    paste0("block", as.integer(block))
  } else {
    block_names
  }

  rows <- list()
  for (bn in use_blocks) {
    if (!bn %in% block_names) next
    blk        <- jackstraw_result[[bn]]
    comp_names <- names(blk)
    use_comps  <- if (!is.null(component)) {
      paste0("comp", as.integer(component))
    } else {
      comp_names
    }
    block_idx <- as.integer(sub("^block", "", bn))

    for (cn in use_comps) {
      if (!cn %in% comp_names) next
      comp_idx <- as.integer(sub("^comp", "", cn))
      cdata    <- blk[[cn]]
      n_feat   <- length(cdata$p_values)

      rows[[length(rows) + 1L]] <- data.frame(
        block         = block_idx,
        component     = comp_idx,
        feature_index = seq_len(n_feat),
        p_value       = cdata$p_values,
        p_adj         = cdata$p_adj,
        significant   = cdata$significant,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    message("No valid block/component combinations found. Returning NULL.")
    return(NULL)
  }
  out <- do.call(rbind, rows)
  out[order(out$block, out$component, out$p_value), ]
}


# ---------------------------------------------------------------------------
# autoplot S3 methods
# ---------------------------------------------------------------------------

#' Automatic ggplot for Rajive objects
#'
#' S3 dispatch for \code{ggplot2::autoplot}.  Delegates to
#' \code{\link{plot_components}} with a sensible default \code{plot_type}.
#'
#' @param object An object of class \code{"rajive"} or
#'   \code{"jackstraw_rajive"}.
#' @param plot_type Character.  Passed to \code{\link{plot_components}}.
#'   For \code{"rajive"} objects the default is \code{"variance"} (requires
#'   \code{blocks} via \code{...}); for \code{"jackstraw_rajive"} objects the
#'   default is \code{"jackstraw_summary"}.
#' @param ... Additional arguments forwarded to \code{\link{plot_components}}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' n <- 40; pks <- c(30, 20)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#' ggplot2::autoplot(fit, plot_type = "pairs")
#'
#' jsr <- jackstraw_rajive(fit, Y$sim_data)
#' ggplot2::autoplot(jsr)
#' }
#'
#' @export
autoplot.rajive <- function(object,
                            plot_type = c("pairs", "density", "top_features",
                                          "component_heatmap", "variance",
                                          "association", "volcano",
                                          "jackstraw_summary", "stability",
                                          "ajive_diagnostic", "rank_threshold",
                                          "bound_distributions"),
                            ...) {
  plot_type <- match.arg(plot_type)
  plot_components(ajive_output = object, plot_type = plot_type, ...)
}

#' @rdname autoplot.rajive
#' @export
autoplot.jackstraw_rajive <- function(object,
                                      plot_type = c("jackstraw_summary",
                                                    "volcano"),
                                      ...) {
  plot_type <- match.arg(plot_type)
  plot_components(jackstraw_result = object, plot_type = plot_type, ...)
}


# ---------------------------------------------------------------------------
# fortify S3 methods
# ---------------------------------------------------------------------------

#' Fortify Rajive objects into tidy data frames
#'
#' S3 dispatch for \code{ggplot2::fortify}.  Returns a tidy
#' \code{data.frame} extracted from the object, suitable for use in
#' \code{ggplot2} pipelines.
#'
#' @param model An object of class \code{"rajive"} or
#'   \code{"jackstraw_rajive"}.
#' @param data Ignored; present for S3 compatibility.
#' @param what Character scalar passed to \code{\link{extract_components}}.
#'   For \code{"rajive"} objects the default is \code{"scores"};
#'   for \code{"jackstraw_rajive"} objects the default is
#'   \code{"significance"}.
#' @param ... Additional arguments forwarded to \code{\link{extract_components}}.
#'
#' @return A \code{data.frame} in long format.
#'
#' @examples
#' \donttest{
#' n <- 40; pks <- c(30, 20)
#' Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
#'                       n = n, pks = pks, dist.type = 1)
#' fit <- Rajive(Y$sim_data, c(4, 3))
#' head(ggplot2::fortify(fit))
#'
#' jsr <- jackstraw_rajive(fit, Y$sim_data)
#' head(ggplot2::fortify(jsr))
#' }
#'
#' @export
fortify.rajive <- function(model, data = NULL,
                           what = c("scores", "loadings", "variance",
                                    "rank_diagnostics"),
                           ...) {
  what <- match.arg(what)
  extract_components(ajive_output = model, what = what, format = "long", ...)
}

#' @rdname fortify.rajive
#' @export
fortify.jackstraw_rajive <- function(model, data = NULL,
                                     what = "significance",
                                     ...) {
  extract_components(what = what, source = "jackstraw",
                     jackstraw_result = model, ...)
}

