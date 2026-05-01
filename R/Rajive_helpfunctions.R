#' Computes the robust SVD of a matrix
#' Using robRsvd
#'
#'
#' @param X Matrix. X matrix.
#' @param rank Integer. Rank of SVD decomposition
#'
#' @return List. The SVD of X.



get_svd_robustH <- function(X, rank=NULL){

  if(is.null(rank)){
    decomposition <- RobRSVD.all(X)
    decomposition
  } else{
    decomposition <- RobRSVD.all(X, nrank = rank)
    decomposition
  }

}


get_sv_threshold <- function(singular_values, rank){

  .5 * (singular_values[rank] + singular_values[rank + 1])
}


#' Truncates a robust SVD.
#'
#' Removes columns from the U, D, V matrix computed form an SVD.
#'
#'
#' @param decomposition List. List with entries 'u', 'd', and 'v'from the svd function.
#' @param rank List. List with entries 'u', 'd', and 'v'from the svd function.
#'
#' @return The trucated robust SVD of X.
truncate_svd <- function(decomposition, rank){

  if(rank==0){
    n <- dim(decomposition[['u']])[1]
    d <- dim(decomposition[['v']])[1]
    decomposition[['u']] <- matrix(0, ncol=1, nrow=n)
    decomposition[['d']] <- 0
    decomposition[['v']] <- matrix(0, ncol=1, nrow=d)
  }else{
    decomposition[['u']] <- decomposition[['u']][, 1:rank, drop=FALSE]
    decomposition[['d']] <- decomposition[['d']][1:rank]
    decomposition[['v']] <- decomposition[['v']][, 1:rank, drop=FALSE]
  }

  decomposition
}

#' Reconstruces the original matrix from its robust SVD.
#'
#' Computes UDV^T to get the approximate (or full) X matrix.
#'
#' @param decomposition List. List with entries 'u', 'd', and 'v'from the svd function.
#'
#' @return Matrix. The original matrix.
svd_reconstruction <- function(decomposition){

  # decomposition rank -- need to truncated singluar values
  r <- dim(decomposition[['u']])[2]

  decomposition[['u']]  %*%
    diag(decomposition[['d']][1:r], nrow=r, ncol=r) %*%
    t(decomposition[['v']])

}

#' Gets the wedin bounds
#'
#' @param X Matrix. The data matrix.
#' @param SVD List. The SVD decomposition of the matrix. List with entries 'u', 'd', and 'v'from the svd function.
#' @param signal_rank Integer.
#' @param num_samples Integer. Number of vectors selected for resampling procedure.


get_wedin_bound_samples <- function(X, SVD, signal_rank, num_samples=1000){

  # resample for U and V
  U_perp <- SVD[['u']][ , -(1:signal_rank)]
  U_sampled_norms <- wedin_bound_resampling(X=X,
                                            perp_basis=U_perp,
                                            right_vectors=FALSE,
                                            num_samples=num_samples)

  V_perp <- SVD[['v']][ , -(1:signal_rank)]
  V_sampled_norms <- wedin_bound_resampling(X=X,
                                            perp_basis=V_perp,
                                            right_vectors=TRUE,
                                            num_samples=num_samples)

  sigma_min <- SVD[['d']][signal_rank]
  wedin_bound_samples <- mapply(function(u, v)  min(max(u, v)/sigma_min, 1)^2, U_sampled_norms, V_sampled_norms)

  wedin_bound_samples
}

#' Resampling procedure for the wedin bound
#'
#' @param X Matrix. The data matrix.
#' @param perp_basis Matrix. Either U_perp or V_perp: the remaining left/right singluar vectors of X after estimating the signal rank.
#' @param right_vectors Boolean. Right multiplication or left multiplication.
#' @param num_samples Integer. Number of vectors selected for resampling procedure.
#' @importFrom foreach %dopar%

wedin_bound_resampling <- function(X, perp_basis, right_vectors, num_samples=1000){

  rank <- dim(perp_basis)[2]
  #resampled_norms <- rep(0, num_samples)
  numCores <- 2
  doParallel::registerDoParallel(numCores)
  resampled_norms <- foreach::foreach (s=1:num_samples) %dopar% {

    sampled_col_index <- sample.int(n=dim(perp_basis)[2],
                                    size=rank,
                                    replace=TRUE)


    perp_resampled <- perp_basis[ , sampled_col_index]

    if(right_vectors){
      resampled_projection <- X %*% perp_resampled
    } else{
      resampled_projection <- t(perp_resampled) %*% X
    }

    # operator L2 norm
    norm(resampled_projection,
         type='2')
  }

  as.numeric(resampled_norms)
}

#' Estimate the wedin bound for a data matrix.
#'
#' Samples from the random direction bound. Returns on the scale of squared singular value.
#'
#' @param n_obs The number of observations.
#' @param dims The number of features in each data matrix
#' @param num_samples Integer. Number of vectors selected for resampling procedure.
#' @importFrom stats rnorm
#'
#' @return rand_dir_samples

get_random_direction_bound_robustH <- function(n_obs, dims, num_samples=1000){

dims1 = as.list(dims)
n_blocks <- length(dims)
numCores <- 2
doParallel::registerDoParallel(numCores)

rand_dir_samples <- foreach::foreach (s=1:num_samples, .export=c("get_svd_robustH", "RobRSVD.all", "RobRSVD1")) %dopar% {

  X <- lapply(dims1, function(l) matrix(rnorm(n_obs * l, mean=0,sd=1), n_obs, l))
  rand_subspaces <- lapply(X, function(l) get_svd_robustH(l)[['u']])

  M <- do.call(cbind, rand_subspaces)
  M_svd <- get_svd_robustH(M, rank=min(dims))

  M_svd[['d']][1]^2

}

as.numeric(rand_dir_samples)
}





#' Block Scores
#'
#' Gets the block scores from the Rajive decomposition
#'
#' @param ajive_output List. The decomposition from Rajive
#' @param k Integer. The index of the data block
#' @param type Character. Joint or individual
#'
#' @return The block scores
#'
#' @examples
#' \donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'get_block_scores(ajive.results.robust, 2, 'joint')
#'}
#' @export


get_block_scores <- function(ajive_output, k, type){
  if(! type  %in% c('joint', 'individual')){
    stop('type must be: joint or individual')
  }
  if (type == 'joint')   res <- ajive_output$block_decomps[[3*(k-1)+2]]$u
  if (type == 'individual')   res <- ajive_output$block_decomps[[3*(k-1)+1]]$u
  
  return(res)

  }


#' Block Loadings
#'
#' Gets the block loadings from the Rajive decomposition
#'
#' @param ajive_output List. The decomposition from Rajive
#' @param k Integer. The index of the data block
#' @param type Character. Joint or individual
#'
#' @return The block loadings
#' @examples 
#'\donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'get_block_loadings(ajive.results.robust, 2, 'joint')
#'}

#'
#' @export
get_block_loadings <- function(ajive_output, k, type){
  if(! type  %in% c('joint', 'individual')){
    stop('type must be: joint or individual')
  }
 if (type == 'joint')   res <- ajive_output$block_decomps[[3*(k-1)+2]]$v
  if (type == 'individual')   res <- ajive_output$block_decomps[[3*(k-1)+1]]$v
return(res)
  }



#' Joint Rank
#'
#' Gets the joint rank from the Rajive decomposition
#'
#' @param ajive_output List. The decomposition from Rajive
#'
#' @return The joint rank
#' @examples  
#'\donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'get_joint_rank(ajive.results.robust)
#'}

#' @export
get_joint_rank <- function(ajive_output){
  ajive_output$joint_rank
}

#' Individual Rank
#'
#' Gets the individual ranks from the Rajive decomposition
#'
#' @param ajive_output List. The decomposition from Rajive
#' @param k Integer. The index of the data block.
#'
#'
#' @return The individual ranks
#'
#'
#' @examples 
#' \donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'get_individual_rank(ajive.results.robust, 2)
#'}

#' @export
get_individual_rank <- function(ajive_output, k){

  ajive_output$block_decomps[[3*(k-1)+1]]$rank
}


#' Decomposition Heatmaps
#'
#' Visualization of the RaJIVE decomposition, it shows heatmaps of the decomposition obtained by RaJIVE
#'
#'
#' @param blocks List. The initial data blocks.
#' @param jive_results_robust List. The RaJIVE decomposition.
#'
#' @return The heatmap of the decomposition
#'
#' @examples
#' \donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'decomposition_heatmaps_robustH(data.ajive, ajive.results.robust)
#'}

#' @export




decomposition_heatmaps_robustH <- function (blocks, jive_results_robust)
{
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("The package 'cowplot' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  K <- length(blocks)


heatmap_listR <- list()

heatmap_listR <- list()
for (k in 1:K) {
  heatmap_listR[[k]] <- data_heatmap(blocks[[k]], ylab = ifelse(k ==
                                                                  1, "observations", ""), show_color_bar = FALSE)
  heatmap_listR[[K + k]] <- data_heatmap(jive_results_robust$block_decomps[[3*(k-1)+2]][["full"]],
                                         ylab = ifelse(k == 1, "joint", ""), show_color_bar = FALSE)
  heatmap_listR[[2 * K + k]] <- data_heatmap(jive_results_robust$block_decomps[[3*(k-1)+1]][["full"]],
                                             ylab = ifelse(k == 1, "individual", ""),
                                             show_color_bar = FALSE)
  heatmap_listR[[3 * K + k]] <- data_heatmap(jive_results_robust$block_decomps[[3*k]],
                                             ylab = ifelse(k == 1, "noise", ""), show_color_bar = FALSE)
}
cowplot::plot_grid(plotlist = heatmap_listR, ncol = K)
}

#' Decomposition Heatmaps
#'
#' Visualization of the RaJIVE decomposition, it shows heatmaps of the decomposition obtained by RaJIVE
#'
#'
#' @param data List. The initial data blocks.
#' @param show_color_bar Boolean.
#' @param title Character.
#' @param xlab Character.
#' @param ylab Character
#'
#' @import ggplot2
#' @importFrom grDevices rainbow


data_heatmap <- function (data, show_color_bar = TRUE, title = "", xlab = "",
                          ylab = "")
{
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("The package 'reshape2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  reshaped_data <- as.data.frame(reshape2::melt(data))
  colnames(reshaped_data) <- c("obs", "var", "value")
  ggplot(data = reshaped_data, aes_string(x = "var",
                                          y = "obs")) +
    geom_raster(aes_string(fill = "value"), show.legend = show_color_bar) +
    scale_fill_gradientn(colours = rainbow(10)) +
    theme(panel.background = element_blank(),  axis.line = element_blank(), legend.position = "bottom") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) + labs(title = title, x = xlab, y = ylab)
}



#' Proportions of variance explained
#'
#' Gets the variance explained by each component of the Rajive decomposition
#'
#' @param ajiveResults List. The decomposition from Rajive
#' @param blocks List. The initial data blocks
#'
#' @return The proportion of variance explained by each component
#' 
#' @examples 
#' \donttest{
#'n <- 10
#'pks <- c(20, 10)
#'Y <- ajive.data.sim(K =2, rankJ = 2, rankA = c(7, 4), n = n,
#'                  pks = pks, dist.type = 1)
#'initial_signal_ranks <-  c(7, 4)
#'data.ajive <- list((Y$sim_data[[1]]), (Y$sim_data[[2]]))
#'ajive.results.robust <- Rajive(data.ajive, initial_signal_ranks)
#'showVarExplained_robust(ajive.results.robust, data.ajive)
#'}
#'
#' @export

showVarExplained_robust <- function(ajiveResults, blocks){
  l <- length(blocks)
  # joint variance
  # joint is the second component for all 3
  VarJoint = rep(0, l)
  for (i in 1:l) VarJoint[i] = norm(as.matrix(ajiveResults$block_decomps[[3*(i-1)+2]][[1]]),
                                    type = "F")^2/norm(blocks[[i]], type = "F")^2

  # individual variances
  # individual is the first component for all 3
  VarIndiv = rep(0, l)
  for (i in 1:l) VarIndiv[i] = norm(as.matrix(ajiveResults$block_decomps[[3*(i-1)+1]][[1]]),
                                    type = "F")^2/norm(blocks[[i]], type = "F")^2

  # residual variance
  VarSubtr = 1 - VarJoint - VarIndiv

  VarProp <- list(VarJoint, VarIndiv, VarSubtr)
  names(VarProp) <- c('Joint', 'Indiv', 'Resid')
  VarProp
}


# ---------------------------------------------------------------------------
# Utility functions for interpreting RaJIVE results
# ---------------------------------------------------------------------------

#' Print method for rajive objects
#'
#' Displays a concise summary of the RaJIVE decomposition: number of blocks,
#' estimated joint rank, and individual rank for each data block.
#'
#' @param x An object of class \code{"rajive"} returned by \code{\link{Rajive}}.
#' @param ... Ignored.
#'
#' @return \code{x} invisibly.
#'
#' @examples
#' \donttest{
#' n <- 30; pks <- c(40, 30)
#' Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = n,
#'                     pks = pks, dist.type = 1)
#' res <- Rajive(Y$sim_data, c(5, 4))
#' print(res)
#' }
#'
#' @export
print.rajive <- function(x, ...) {
  K <- length(x$block_decomps) / 3L
  cat("RaJIVE Decomposition\n")
  cat(sprintf("  Number of blocks : %d\n", K))
  cat(sprintf("  Joint rank       : %d\n", x$joint_rank))
  indiv_ranks <- vapply(seq_len(K),
                        function(k) get_individual_rank(x, k),
                        numeric(1L))
  cat(sprintf("  Individual ranks : %s\n",
              paste(indiv_ranks, collapse = ", ")))
  invisible(x)
}


#' Summary method for rajive objects
#'
#' Returns (and prints) a \code{data.frame} with the joint rank and
#' individual rank for every data block.
#'
#' @param object An object of class \code{"rajive"} returned by
#'   \code{\link{Rajive}}.
#' @param ... Ignored.
#'
#' @return A \code{data.frame} with columns \code{block},
#'   \code{joint_rank}, and \code{individual_rank}, returned invisibly.
#'
#' @examples
#' \donttest{
#' n <- 30; pks <- c(40, 30)
#' Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = n,
#'                     pks = pks, dist.type = 1)
#' res <- Rajive(Y$sim_data, c(5, 4))
#' summary(res)
#' }
#'
#' @export
summary.rajive <- function(object, ...) {
  df <- get_all_ranks(object)
  print(df, row.names = FALSE)
  invisible(df)
}


#' Joint Scores
#'
#' Returns the shared (cross-block) joint score matrix from a RaJIVE
#' decomposition.  Each column corresponds to one joint component and each
#' row to one observation.
#'
#' @param ajive_output List returned by \code{\link{Rajive}}.
#'
#' @return An \eqn{n \times \text{joint\_rank}} numeric matrix of joint scores.
#'
#' @examples
#' \donttest{
#' n <- 30; pks <- c(40, 30)
#' Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = n,
#'                     pks = pks, dist.type = 1)
#' res <- Rajive(Y$sim_data, c(5, 4))
#' get_joint_scores(res)
#' }
#'
#' @export
get_joint_scores <- function(ajive_output) {
  ajive_output[["joint_scores"]]
}


#' Extract a reconstructed block matrix
#'
#' Returns the full reconstructed matrix for the joint (\eqn{J}), individual
#' (\eqn{I}), or noise (\eqn{E}) component of a single data block from a
#' RaJIVE decomposition.
#'
#' @param ajive_output List returned by \code{\link{Rajive}}.
#' @param k Positive integer; index of the data block.
#' @param type Character string; one of \code{"joint"}, \code{"individual"},
#'   or \code{"noise"}.
#'
#' @return The reconstructed matrix for the requested component and block.
#'   Returns \code{NA} if \code{\link{Rajive}} was called with
#'   \code{full = FALSE}.
#'
#' @examples
#' \donttest{
#' n <- 30; pks <- c(40, 30)
#' Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = n,
#'                     pks = pks, dist.type = 1)
#' res <- Rajive(Y$sim_data, c(5, 4))
#' J1 <- get_block_matrix(res, k = 1, type = "joint")
#' }
#'
#' @export
get_block_matrix <- function(ajive_output, k, type = c("joint", "individual", "noise")) {
  type <- match.arg(type)
  if (type == "joint")       return(ajive_output$block_decomps[[3L * (k - 1L) + 2L]][["full"]])
  if (type == "individual")  return(ajive_output$block_decomps[[3L * (k - 1L) + 1L]][["full"]])
  if (type == "noise")       return(ajive_output$block_decomps[[3L * k]])
}


#' Summary table of all ranks
#'
#' Returns a \code{data.frame} with the joint rank and individual rank for
#' every data block, making it easy to inspect all estimated ranks at once.
#'
#' @param ajive_output List returned by \code{\link{Rajive}}.
#'
#' @return A \code{data.frame} with columns \code{block},
#'   \code{joint_rank}, and \code{individual_rank}.
#'
#' @examples
#' \donttest{
#' n <- 30; pks <- c(40, 30)
#' Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = n,
#'                     pks = pks, dist.type = 1)
#' res <- Rajive(Y$sim_data, c(5, 4))
#' get_all_ranks(res)
#' }
#'
#' @export
get_all_ranks <- function(ajive_output) {
  K <- length(ajive_output$block_decomps) / 3L
  indiv_ranks <- vapply(seq_len(K),
                        function(k) get_individual_rank(ajive_output, k),
                        numeric(1L))
  data.frame(
    block           = paste0("block", seq_len(K)),
    joint_rank      = rep(ajive_output$joint_rank, K),
    individual_rank = indiv_ranks,
    stringsAsFactors = FALSE
  )
}


#' Bar chart of variance explained
#'
#' Produces a stacked bar chart showing the proportion of total variance
#' explained by the joint, individual, and residual components for each
#' data block.
#'
#' @param ajive_output List returned by \code{\link{Rajive}}.
#' @param blocks List of data matrices (the same list passed to
#'   \code{\link{Rajive}}).
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \donttest{
#' n <- 30; pks <- c(40, 30)
#' Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(5, 4), n = n,
#'                     pks = pks, dist.type = 1)
#' data.ajive <- Y$sim_data
#' res <- Rajive(data.ajive, c(5, 4))
#' plot_variance_explained(res, data.ajive)
#' }
#'
#' @import ggplot2
#' @export
plot_variance_explained <- function(ajive_output, blocks) {
  var_exp <- showVarExplained_robust(ajive_output, blocks)
  K <- length(blocks)

  plot_df <- data.frame(
    block     = rep(paste0("block", seq_len(K)), 3L),
    component = c(rep("Joint",      K),
                  rep("Individual", K),
                  rep("Residual",   K)),
    proportion = c(unlist(var_exp[["Joint"]]),
                   unlist(var_exp[["Indiv"]]),
                   unlist(var_exp[["Resid"]])),
    stringsAsFactors = FALSE
  )
  plot_df$component <- factor(plot_df$component,
                               levels = c("Residual", "Individual", "Joint"))

  ggplot(plot_df, aes(x = block, y = proportion, fill = component)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(
      "Joint"      = "#2c7bb6",
      "Individual" = "#abd9e9",
      "Residual"   = "#d7191c"
    )) +
    labs(
      title = "Proportion of Variance Explained",
      x     = "Data Block",
      y     = "Proportion",
      fill  = "Component"
    ) +
    ylim(0, 1) +
    theme_bw()
}


#' Scatter plot of block scores
#'
#' Plots two score components against each other for a given data block and
#' component type (joint or individual), making it easy to visualise the
#' latent structure captured by the RaJIVE decomposition.
#'
#' @param ajive_output List returned by \code{\link{Rajive}}.
#' @param k Positive integer; index of the data block.
#' @param type Character string; \code{"joint"} or \code{"individual"}.
#' @param comp_x Positive integer; index of the component to plot on the
#'   x-axis.  Default \code{1}.
#' @param comp_y Positive integer; index of the component to plot on the
#'   y-axis.  Default \code{2}.
#' @param group Optional factor or vector (length equal to the number of
#'   observations) used to colour the points.  \code{NULL} (default) gives
#'   uniform colouring.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \donttest{
#' n <- 50; pks <- c(60, 40)
#' Y <- ajive.data.sim(K = 2, rankJ = 3, rankA = c(5, 4), n = n,
#'                     pks = pks, dist.type = 1)
#' res <- Rajive(Y$sim_data, c(5, 4))
#' plot_scores(res, k = 1, type = "joint")
#' plot_scores(res, k = 2, type = "individual", comp_x = 1, comp_y = 2)
#' }
#'
#' @import ggplot2
#' @export
plot_scores <- function(ajive_output, k,
                        type   = c("joint", "individual"),
                        comp_x = 1L,
                        comp_y = 2L,
                        group  = NULL) {
  type   <- match.arg(type)
  scores <- get_block_scores(ajive_output, k, type)

  n_comp <- ncol(scores)
  if (comp_x > n_comp || comp_y > n_comp) {
    stop(sprintf(
      "'comp_x' and 'comp_y' must each be <= %d (the number of %s components for block %d).",
      n_comp, type, k
    ), call. = FALSE)
  }

  plot_df <- data.frame(x = scores[, comp_x], y = scores[, comp_y])

  if (!is.null(group)) {
    if (length(group) != nrow(scores)) {
      stop("'group' must have the same length as the number of observations.",
           call. = FALSE)
    }
    plot_df$group <- as.factor(group)
  }

  type_label <- paste0(toupper(substring(type, 1L, 1L)), substring(type, 2L))
  p <- ggplot(plot_df, aes(x = x, y = y))

  if (!is.null(group)) {
    p <- p + geom_point(aes(colour = group))
  } else {
    p <- p + geom_point()
  }

  p +
    labs(
      title = sprintf("%s scores: block %d (comp %d vs comp %d)",
                      type_label, k, comp_x, comp_y),
      x = sprintf("Component %d", comp_x),
      y = sprintf("Component %d", comp_y)
    ) +
    theme_bw()
}
