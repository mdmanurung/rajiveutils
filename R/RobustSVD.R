

#' Computes the robust SVD of a matrix
#'
#' Thin R wrapper around the RcppArmadillo implementation \code{RobRSVD_all_cpp}.
#' The interface and return value are identical to the original pure-R version.
#'
#' @param data Matrix. X matrix.
#' @param nrank Integer. Rank of SVD decomposition
#' @param svdinit List or NULL. Optional classical SVD used for initialization.
#'   If NULL (default), computed lazily after input validation.
#' @param weights Optional numeric or logical matrix with the same dimensions
#'   as \code{data}. Positive / \code{TRUE} entries are treated as observed;
#'   zero / \code{FALSE} entries are ignored by the weighted SVD path.
#' @importFrom stats median
#' @return List with entries \code{d}, \code{u}, \code{v}.  When
#'   \code{nrank <= 0}, returns \code{numeric(0)} singular values and
#'   zero-column \code{u}/\code{v} matrices without calling the C++ backend.

RobRSVD.all <- function(data, nrank = min(dim(data)), svdinit = NULL,
                        weights = NULL)
{
  if (!is.matrix(data)) {
    cli::cli_abort(
      c("`data` must be a numeric matrix."),
      class = "rajiveplus_invalid_input"
    )
  }
  weights <- .normalize_robsvd_weights(data, weights)
  if (is.null(weights) && !all(is.finite(data))) {
    cli::cli_abort(
      c("`data` contains non-finite values (NA/NaN/Inf)."),
      class = "rajiveplus_invalid_input"
    )
  }
  if (!is.null(weights) && any(weights & !is.finite(data))) {
    cli::cli_abort(
      c("Observed entries of `data` must be finite when `weights` are supplied."),
      class = "rajiveplus_invalid_input"
    )
  }
  nrank <- as.integer(nrank)
  if (length(nrank) != 1L || is.na(nrank)) {
    cli::cli_abort(
      c("`nrank` must be a single integer."),
      class = "rajiveplus_invalid_input"
    )
  }
  if (nrank <= 0L) {
    return(list(
      d = numeric(0),
      u = matrix(0, nrow = nrow(data), ncol = 0L),
      v = matrix(0, nrow = ncol(data), ncol = 0L)
    ))
  }
  if (!is.null(weights)) {
    if (all(weights)) {
      if (!all(is.finite(data))) {
        cli::cli_abort(
          c("`data` contains non-finite observed values."),
          class = "rajiveplus_invalid_input"
        )
      }
    } else {
      return(.RobRSVD_all_weighted_R(data, weights, nrank = nrank))
    }
  }
  if (is.null(svdinit)) {
    svdinit <- svd(data)
  }

  RobRSVD_all_cpp(
    data   = data,
    nrank  = nrank,
    sinit1 = svdinit$d[1],
    uinit1 = svdinit$u[, 1, drop = TRUE],
    vinit1 = svdinit$v[, 1, drop = TRUE]
  )
}

.normalize_robsvd_weights <- function(data, weights) {
  if (is.null(weights)) return(NULL)
  if (!is.matrix(weights) || !identical(dim(weights), dim(data))) {
    cli::cli_abort(
      "`weights` must be a matrix with the same dimensions as `data`.",
      class = "rajiveplus_invalid_input"
    )
  }
  if (!is.numeric(weights) && !is.logical(weights)) {
    cli::cli_abort(
      "`weights` must be numeric or logical.",
      class = "rajiveplus_invalid_input"
    )
  }
  if (anyNA(weights) || any(!is.finite(as.numeric(weights)))) {
    cli::cli_abort(
      "`weights` must contain only finite, non-missing values.",
      class = "rajiveplus_invalid_input"
    )
  }
  weights > 0
}

.RobRSVD_all_weighted_R <- function(data, weights, nrank, max_iter = 100L,
                                    tol = 1e-7) {
  r <- min(as.integer(nrank), min(dim(data)))
  if (r <= 0L || !any(weights)) {
    return(list(
      d = numeric(0),
      u = matrix(0, nrow = nrow(data), ncol = 0L),
      v = matrix(0, nrow = ncol(data), ncol = 0L)
    ))
  }

  observed <- weights
  safe <- data
  safe[!observed | !is.finite(safe)] <- NA_real_
  filled <- safe
  for (j in seq_len(ncol(filled))) {
    replacement <- mean(filled[observed[, j], j], na.rm = TRUE)
    if (!is.finite(replacement)) replacement <- 0
    filled[!observed[, j] | !is.finite(filled[, j]), j] <- replacement
  }
  filled[!is.finite(filled)] <- 0

  previous_objective <- Inf
  for (iter in seq_len(max_iter)) {
    sv <- svd(filled, nu = r, nv = r)
    u <- sv$u[, seq_len(r), drop = FALSE]
    d <- sv$d[seq_len(r)]
    v <- sv$v[, seq_len(r), drop = FALSE]
    recon <- u %*% (diag(d, nrow = r, ncol = r) %*% t(v))
    objective <- sum((data[observed] - recon[observed])^2)
    filled[!observed] <- recon[!observed]
    filled[observed] <- data[observed]
    if (is.finite(previous_objective) &&
        abs(previous_objective - objective) <= tol * (abs(previous_objective) + tol)) {
      break
    }
    previous_objective <- objective
  }

  sv <- svd(filled, nu = r, nv = r)
  u <- sv$u[, seq_len(r), drop = FALSE]
  d <- sv$d[seq_len(r)]
  v <- sv$v[, seq_len(r), drop = FALSE]

  fully_masked_rows <- rowSums(observed) == 0L
  fully_masked_cols <- colSums(observed) == 0L
  if (any(fully_masked_rows)) u[fully_masked_rows, ] <- 0
  if (any(fully_masked_cols)) v[fully_masked_cols, ] <- 0

  list(d = d, u = u, v = v)
}


#' Single robust rank-1 component (Rcpp wrapper)
#'
#' Thin R wrapper around \code{RobRSVD1_cpp}, preserved for backward
#' compatibility and direct testing.
#'
#' @param data Matrix. X matrix.
#' @param huberk Numeric. Huber k tuning constant.
#' @param niter Integer. Maximum iterations.
#' @param tol Numeric. Convergence tolerance.
#' @param sinit Numeric. Initial singular value.
#' @param uinit Numeric vector. Initial left singular vector.
#' @param vinit Numeric vector. Initial right singular vector.
#' @return List with entries \code{s}, \code{u}, \code{v}.

RobRSVD1 <- function(data, huberk = 1.345, niter = 1000,
                     tol = 1e-05, sinit, uinit, vinit)
{
  RobRSVD1_cpp(
    data   = data,
    sinit  = sinit,
    uinit  = uinit,
    vinit  = vinit,
    huberk = huberk,
    niter  = as.integer(niter),
    tol    = tol
  )
}
