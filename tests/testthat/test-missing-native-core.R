test_that("incomplete core preserves complete-data fixed-rank decomposition", {
  set.seed(8201)
  Y <- ajive.data.sim(K = 2, rankJ = 1, rankA = c(1, 1),
                      n = 12, pks = c(6, 5), dist.type = 1)
  blocks <- Y$sim_data
  ranks <- c(2L, 2L)

  complete <- rajiveplus:::.Rajive_core(
    blocks, ranks, joint_rank = 1L,
    n_wedin_samples = NA, n_rand_dir_samples = NA,
    num_cores = 1L, identifiability_norm = "l2",
    rank_only = FALSE
  )
  native <- rajiveplus:::.Rajive_incomplete(
    blocks, ranks, joint_rank = 1L,
    mask = lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x))),
    missing_control = rajive_missing_control(center = FALSE, scale = FALSE,
                                             normalize = FALSE)
  )

  expect_equal(native$joint_rank, complete$joint_rank)
  expect_equal(abs(crossprod(native$joint_scores, complete$joint_scores)),
               diag(1), tolerance = 1e-6)
  expect_equal(get_block_matrix(native, 1, "joint"),
               get_block_matrix(complete, 1, "joint"),
               tolerance = 1e-6)
})

test_that("incomplete core is invariant to masked numeric values", {
  set.seed(8202)
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[2, 3] <- FALSE
  mask$block2[5, ] <- FALSE

  corrupt <- blocks
  corrupt$block1[2, 3] <- 1e7
  corrupt$block2[5, ] <- -1e7

  a <- rajiveplus:::.Rajive_incomplete(blocks, c(2L, 2L), joint_rank = 1L,
                                       mask = mask)
  b <- rajiveplus:::.Rajive_incomplete(corrupt, c(2L, 2L), joint_rank = 1L,
                                       mask = mask)

  expect_equal(abs(crossprod(a$joint_scores, b$joint_scores)),
               diag(1), tolerance = 1e-6)
  expect_equal(get_block_matrix(a, 1, "joint")[mask$block1],
               get_block_matrix(b, 1, "joint")[mask$block1],
               tolerance = 1e-6)
  expect_equal(get_block_matrix(a, 2, "individual")[mask$block2],
               get_block_matrix(b, 2, "individual")[mask$block2],
               tolerance = 1e-6)
})

test_that("incomplete core recovers a simple joint subspace under sample-block missingness", {
  set.seed(8203)
  n <- 18
  z <- scale(seq_len(n), center = TRUE, scale = TRUE)
  X1 <- z %*% matrix(c(2, -1, 0.5, 1), nrow = 1) +
    matrix(rnorm(n * 4, sd = 0.05), n)
  X2 <- z %*% matrix(c(-1, 1.5, 0.25), nrow = 1) +
    matrix(rnorm(n * 3, sd = 0.05), n)
  blocks <- list(block1 = X1, block2 = X2)
  complete <- Rajive(blocks, c(1L, 1L), joint_rank = 1L,
                     n_wedin_samples = NA, n_rand_dir_samples = NA)

  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block2[seq(2, n, by = 4), ] <- FALSE
  incomplete <- rajiveplus:::.Rajive_incomplete(blocks, c(1L, 1L),
                                                joint_rank = 1L, mask = mask)

  alignment <- abs(cor(complete$joint_scores[, 1],
                       incomplete$joint_scores[, 1]))
  expect_gt(alignment, 0.8)
})

test_that("incomplete core records finite convergence diagnostics", {
  set.seed(8204)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[1, 1] <- FALSE

  fit <- rajiveplus:::.Rajive_incomplete(blocks, c(2L, 2L),
                                         joint_rank = 1L, mask = mask)

  expect_true(is.finite(fit$missing$convergence$objective))
  expect_equal(length(fit$missing$convergence$block_objective), 2L)
})

test_that("native missing fits support full false with missing cells", {
  set.seed(8205)
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[2, 3] <- FALSE
  mask$block2[5, ] <- FALSE

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L, full = FALSE)
  diag <- get_missing_diagnostics(fit)
  recon <- get_reconstructed_blocks(fit, type = "joint_individual",
                                    scale = "standardized")

  expect_s3_class(fit, "rajive_incomplete")
  expect_true(is.na(get_block_matrix(fit, 1, "joint")))
  expect_true(all(is.finite(diag$variance_explained$Joint)))
  expect_equal(dim(recon$block1), dim(blocks$block1))
  expect_true(all(is.finite(recon$block1)))
})

test_that("all-observed native branch supports full false", {
  set.seed(8206)
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))

  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L, full = FALSE,
                n_wedin_samples = NA, n_rand_dir_samples = NA)
  diag <- get_missing_diagnostics(fit)
  recon <- get_reconstructed_blocks(fit, type = "joint_individual")

  expect_s3_class(fit, "rajive_incomplete")
  expect_true(is.na(get_block_matrix(fit, 1, "joint")))
  expect_true(all(is.finite(diag$variance_explained$Joint)))
  expect_equal(dim(recon$block2), dim(blocks$block2))
  expect_true(all(is.finite(recon$block2)))
})
