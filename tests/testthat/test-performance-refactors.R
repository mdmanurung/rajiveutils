library(testthat)
library(rajiveplus)

old_ols_f_stat_matrix_perf <- function(Y_t, x) {
  n <- length(x)
  x_c <- x - mean(x)
  row_means <- rowMeans(Y_t)
  Y_c <- Y_t - row_means
  xTx <- sum(x_c^2)
  beta1 <- as.numeric(Y_c %*% x_c) / xTx
  E <- Y_c - tcrossprod(beta1, x_c)
  sse1 <- rowSums(E^2)
  ss0 <- rowSums(Y_c^2)
  F_stat <- (ss0 - sse1) / (sse1 / (n - 2))
  is_constant <- (rajiveplus:::rowVars_fast(Y_t) < .Machine$double.eps)
  F_stat[is_constant] <- NA_real_
  F_stat
}

old_generate_null_f_stats_perf <- function(X_t, joint_comp_scores, n_null) {
  d <- nrow(X_t)
  sampled_idx <- matrix(
    sample.int(d, size = d * n_null, replace = TRUE),
    nrow = d, ncol = n_null
  )
  null_f <- matrix(NA_real_, nrow = d, ncol = n_null)
  for (s in seq_len(n_null)) {
    rows <- X_t[sampled_idx[, s], , drop = FALSE]
    scores_perm <- sample(joint_comp_scores)
    null_f[, s] <- old_ols_f_stat_matrix_perf(rows, scores_perm)
  }
  null_f
}

leading_sv_svd_perf <- function(x) {
  d <- svd(x, nu = 0L, nv = 0L)$d
  if (length(d) == 0L) 0 else d[[1L]]
}

leading_sv_crossprod_perf <- function(x) {
  if (length(x) == 0L) return(0)
  gram <- if (nrow(x) >= ncol(x)) crossprod(x) else tcrossprod(x)
  vals <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
  sqrt(max(vals[[1L]], 0))
}

time_min_perf <- function(expr, reps = 3L) {
  expr <- substitute(expr)
  env <- parent.frame()
  timings <- numeric(reps)
  for (i in seq_len(reps)) {
    invisible(gc())
    timings[[i]] <- system.time(eval(expr, env))[["elapsed"]]
  }
  min(timings)
}

test_that("PERF-001 leading singular value crossproduct prototype matches SVD", {
  set.seed(1401)
  mats <- list(
    matrix(rnorm(80 * 5), 80, 5),
    matrix(rnorm(5 * 80), 5, 80),
    matrix(rnorm(30 * 30), 30, 30),
    matrix(rnorm(12), 12, 1)
  )
  for (x in mats) {
    expect_equal(
      leading_sv_crossprod_perf(x),
      leading_sv_svd_perf(x),
      tolerance = 1e-10
    )
  }
})

test_that("PERF-002 ols_f_stat_matrix matches frozen pre-refactor formula", {
  set.seed(1402)
  X_t <- matrix(rnorm(120 * 35), 120, 35)
  X_t[1L, ] <- 1
  scores <- rnorm(35)

  expect_equal(
    rajiveplus:::ols_f_stat_matrix(X_t, scores),
    old_ols_f_stat_matrix_perf(X_t, scores),
    tolerance = 1e-10
  )
})

test_that("PERF-002 null F-stat generation preserves frozen RNG/result contract", {
  set.seed(1403)
  X_t <- matrix(rnorm(160 * 40), 160, 40)
  X_t[2L, ] <- 3
  scores <- rnorm(40)

  set.seed(1404)
  old <- old_generate_null_f_stats_perf(X_t, scores, n_null = 8L)
  set.seed(1404)
  cur <- rajiveplus:::generate_null_f_stats(X_t, scores, n_null = 8L)

  expect_equal(cur, old, tolerance = 1e-10)
})

test_that("PERF-002 null F-stat generation is faster than frozen implementation", {
  skip_if_not(
    identical(Sys.getenv("RAJIVE_RUN_PERF"), "1"),
    "set RAJIVE_RUN_PERF=1 to run performance gates"
  )

  set.seed(1405)
  X_t <- matrix(rnorm(2500 * 45), 2500, 45)
  scores <- rnorm(45)
  n_null <- 20L

  set.seed(1406)
  old <- old_generate_null_f_stats_perf(X_t, scores, n_null)
  set.seed(1406)
  cur <- rajiveplus:::generate_null_f_stats(X_t, scores, n_null)
  expect_equal(cur, old, tolerance = 1e-10)

  invisible(rajiveplus:::generate_null_f_stats(X_t, scores, 2L))
  invisible(old_generate_null_f_stats_perf(X_t, scores, 2L))
  set.seed(1407)
  old_time <- time_min_perf(
    old_generate_null_f_stats_perf(X_t, scores, n_null),
    reps = 3L
  )
  set.seed(1407)
  cur_time <- time_min_perf(
    rajiveplus:::generate_null_f_stats(X_t, scores, n_null),
    reps = 3L
  )

  max_ratio <- as.numeric(Sys.getenv("RAJIVE_PERF_MAX_RATIO", "0.85"))
  expect_lt(cur_time / old_time, max_ratio)
})

test_that("PERF-003 fixed-rank bootstrap is treated as opt-in non-equivalent strategy", {
  set.seed(1408)
  Y <- ajive.data.sim(K = 2, rankJ = 1, rankA = c(3, 3), n = 35,
                      pks = c(18, 16), dist.type = 1)
  fit <- Rajive(Y$sim_data, c(3, 3), n_wedin_samples = 10,
                n_rand_dir_samples = 10)

  set.seed(1409)
  fixed_a <- rajiveplus:::.rajive_bootstrap(
    fit, Y$sim_data, c(3, 3), B = 3L, keep = "loadings",
    n_wedin_samples = 10, n_rand_dir_samples = 10,
    joint_rank = fit$joint_rank
  )
  set.seed(1409)
  fixed_b <- rajiveplus:::.rajive_bootstrap(
    fit, Y$sim_data, c(3, 3), B = 3L, keep = "loadings",
    n_wedin_samples = 10, n_rand_dir_samples = 10,
    joint_rank = fit$joint_rank
  )

  expect_equal(fixed_a, fixed_b, tolerance = 1e-10)
})

test_that("PERF-004 preallocated CI row prototype matches list-rbind output", {
  block <- rep(c("block1", "block2"), each = 4L)
  component <- rep(rep(1:2, each = 2L), 2L)
  feature <- paste0("f", seq_along(block))
  estimate <- seq_along(block) / 10
  lower <- estimate - 0.01
  upper <- estimate + 0.01

  old <- do.call(rbind, lapply(seq_along(block), function(i) {
    data.frame(
      target = "loadings", block = block[[i]], component = component[[i]],
      feature = feature[[i]], sample = NA_character_,
      estimate = estimate[[i]], lower = lower[[i]], upper = upper[[i]],
      level = 0.95, method = "percentile", n_replicates = 25L,
      stringsAsFactors = FALSE
    )
  }))
  rownames(old) <- NULL

  cur <- data.frame(
    target = rep("loadings", length(block)),
    block = block,
    component = component,
    feature = feature,
    sample = rep(NA_character_, length(block)),
    estimate = estimate,
    lower = lower,
    upper = upper,
    level = rep(0.95, length(block)),
    method = rep("percentile", length(block)),
    n_replicates = rep(25L, length(block)),
    stringsAsFactors = FALSE
  )

  expect_equal(cur, old, tolerance = 0)
})
