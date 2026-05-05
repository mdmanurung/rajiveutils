# visualization.R — Interpretation and diagnostic helpers for rajiveplus
#
# Public API:
#   extract_components()   — extract tidy data / diagnostic payload from Rajive output
#   plot_components()      — unified plotting entry point
#
# Internal helpers (not exported):
#   .extract_rank_diagnostics()
#   .plot_rank_threshold()
#   .plot_bound_distributions()
#   .plot_ajive_diagnostic()


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
#' @param what Character scalar. What to extract.  Only
#'   \code{"rank_diagnostics"} is currently implemented.
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
extract_components <- function(ajive_output,
                               what   = "rank_diagnostics",
                               format = c("wide", "long"),
                               include_meta = FALSE,
                               ...) {
  format <- match.arg(format)

  if (!identical(what, "rank_diagnostics")) {
    stop(sprintf(
      "what = '%s' is not supported. Currently only 'rank_diagnostics' is implemented.",
      what), call. = FALSE)
  }

  if (!inherits(ajive_output, "rajive")) {
    stop("ajive_output must be of class 'rajive'.", call. = FALSE)
  }

  .extract_rank_diagnostics(ajive_output, format = format, include_meta = include_meta)
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
#'     classification.}
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
#' @param plot_type Character.  One of \code{"ajive_diagnostic"},
#'   \code{"rank_threshold"}, or \code{"bound_distributions"}.
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
                            plot_type = c("ajive_diagnostic",
                                          "rank_threshold",
                                          "bound_distributions"),
                            ...) {
  plot_type <- match.arg(plot_type)
  dots <- list(...)

  # Resolve diagnostics payload
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

  switch(plot_type,
    ajive_diagnostic    = .plot_ajive_diagnostic(diag, dots),
    rank_threshold      = .plot_rank_threshold(diag, dots),
    bound_distributions = .plot_bound_distributions(diag, dots)
  )
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

  if (!has_wedin && !has_random) {
    warning(paste0(
      "Both Wedin and random-direction bounds are absent from joint_rank_sel. ",
      "Diagnostic visualizations will show observed singular values only."),
      call. = FALSE)
  }

  overall_sv_sq_threshold <- jrs[["overall_sv_sq_threshold"]]
  if (is.null(overall_sv_sq_threshold)) overall_sv_sq_threshold <- NA_real_

  cutoff_rule <- if (has_wedin && has_random) {
    "max(wedin, random)"
  } else if (has_wedin) {
    "wedin_only"
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
    wedin_cutoff            = wedin_cutoff,
    rand_cutoff             = rand_cutoff,
    wedin_percentile        = wedin_percentile,
    rand_percentile         = rand_percentile,
    identif_dropped         = identif_dropped,
    cutoff_rule             = cutoff_rule,
    has_wedin               = has_wedin,
    has_random              = has_random
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
    stringsAsFactors        = FALSE
  )
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
  threshold  <- diag$overall_sv_sq_threshold

  subtitle <- if (has_wedin && has_random) {
    sprintf("Joint rank: %d | Threshold: %.4f (%s)",
            diag$joint_rank_estimate, threshold, diag$cutoff_rule)
  } else if (has_wedin) {
    sprintf("Joint rank: %d | Wedin bound only (rand_dir bound unavailable)",
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

  # Both distributions: combine side-by-side
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(paste0(
      "Package 'patchwork' is required to combine both bound-distribution plots. ",
      "Install with: install.packages('patchwork')."),
      call. = FALSE)
  }
  plots[["wedin"]] + plots[["rand_dir"]]
}


# ---------------------------------------------------------------------------
# Internal: .plot_ajive_diagnostic
# ---------------------------------------------------------------------------

.plot_ajive_diagnostic <- function(diag, dots) {
  p_thresh <- .plot_rank_threshold(diag, dots)

  has_wedin  <- diag$has_wedin
  has_random <- diag$has_random

  # No bounds at all: return single panel with annotation already in subtitle
  if (!has_wedin && !has_random) {
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

  bottom <- if (length(dist_plots) == 2L) {
    dist_plots[["wedin"]] + dist_plots[["rand_dir"]]
  } else {
    dist_plots[[1L]]
  }

  p_thresh / bottom
}
