

#' Computes the robust SVD of a matrix
#'
#' Thin R wrapper around the RcppArmadillo implementation \code{RobRSVD_all_cpp}.
#' The interface and return value are identical to the original pure-R version.
#'
#' @param data Matrix. X matrix.
#' @param nrank Integer. Rank of SVD decomposition
#' @param svdinit List or NULL. Optional classical SVD used for initialization.
#'   If NULL (default), computed lazily after input validation.
#' @importFrom stats median
#' @return List with entries \code{d}, \code{u}, \code{v}.  When
#'   \code{nrank <= 0}, returns \code{numeric(0)} singular values and
#'   zero-column \code{u}/\code{v} matrices without calling the C++ backend.

RobRSVD.all <- function(data, nrank = min(dim(data)), svdinit = NULL)
{
  if (!is.matrix(data)) {
    cli::cli_abort(
      c("`data` must be a numeric matrix."),
      class = "rajiveplus_invalid_input"
    )
  }
  if (!all(is.finite(data))) {
    cli::cli_abort(
      c("`data` contains non-finite values (NA/NaN/Inf)."),
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
