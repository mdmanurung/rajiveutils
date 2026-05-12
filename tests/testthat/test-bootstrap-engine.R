test_that(".bootstrap_resample_indices handles observation bootstrap", {
  set.seed(1)
  idx <- rajiveplus:::.bootstrap_resample_indices(n = 10L, sample_frac = 0.6)

  expect_length(idx, 6L)
  expect_true(all(idx >= 1L & idx <= 10L))
})

test_that(".bootstrap_resample_indices handles whole-cluster bootstrap", {
  set.seed(2)
  cluster <- rep(letters[1:4], each = 3L)
  idx <- rajiveplus:::.bootstrap_resample_indices(
    n = length(cluster),
    sample_frac = 1,
    cluster = cluster,
    resample = "cluster"
  )

  counts <- table(cluster[idx])
  expect_true(all(as.integer(counts) %% 3L == 0L))
  expect_equal(length(idx), length(cluster))
})

test_that(".bootstrap_resample_indices respects cluster strata", {
  set.seed(3)
  cluster <- rep(letters[1:6], each = 2L)
  strata <- rep(c("case", "control"), each = 6L)
  idx <- rajiveplus:::.bootstrap_resample_indices(
    n = length(cluster),
    sample_frac = 1,
    cluster = cluster,
    strata = strata,
    resample = "cluster"
  )

  expect_equal(unname(table(strata[idx])), unname(table(strata)))
})

test_that("bootstrap score helpers extract and validate selected score families", {
  fx <- make_small_rajive_fixture()

  expect_equal(
    rajiveplus:::.extract_score_matrix(fx$fit, "joint", NULL),
    fx$fit$joint_scores
  )
  expect_equal(
    rajiveplus:::.extract_score_matrix(fx$fit, "individual", 1L),
    fx$fit$block_decomps[[1L]]$u
  )
  expect_error(
    rajiveplus:::.validate_bootstrap_score_request("individual", NULL, 2L),
    "`score_block`"
  )
  expect_error(
    rajiveplus:::.validate_bootstrap_score_request("individual", 3L, 2L),
    "`score_block`"
  )
})

test_that(".align_and_scatter_scores handles rotations and duplicate rows", {
  ref <- cbind(seq_len(6), c(2, -1, 1, 0, 3, -2))
  Q <- matrix(c(0, 1, -1, 0), nrow = 2)
  idx <- c(1L, 2L, 2L, 3L, 4L, 5L)
  bs <- ref[idx, , drop = FALSE] %*% t(Q)

  aligned <- rajiveplus:::.align_and_scatter_scores(ref, bs, idx, n_ref = 6L)

  expect_equal(dim(aligned$scores), c(6L, 2L))
  expect_true(all(is.finite(aligned$scores[1:5, ])))
  expect_true(all(is.na(aligned$scores[6L, ])))
  expect_gt(abs(stats::cor(ref[1:5, 1L], aligned$scores[1:5, 1L])), 0.99)
  expect_gt(abs(stats::cor(ref[1:5, 2L], aligned$scores[1:5, 2L])), 0.99)

  too_few <- rajiveplus:::.align_and_scatter_scores(
    ref, bs[1:2, , drop = FALSE], idx[1:2], n_ref = 6L
  )
  expect_null(too_few)
})

test_that(".rajive_bootstrap returns requested replicate payloads", {
  fx <- make_small_rajive_fixture()

  set.seed(4)
  reps <- rajiveplus:::.rajive_bootstrap(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    B = 2L,
    sample_frac = 0.75,
    keep = c("loadings", "scores", "joint_rank", "component_cors",
             "indices", "var_explained"),
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  expect_named(reps, c("loadings", "scores", "joint_rank", "component_cors",
                       "indices", "var_explained"))
  expect_length(reps$loadings, 2L)
  expect_equal(dim(reps$loadings[[1L]]), c(8L, 1L, 2L))
  expect_equal(dim(reps$scores), c(12L, 1L, 2L))
  expect_length(reps$joint_rank, 2L)
  expect_equal(dim(reps$component_cors), c(2L, 1L))
  expect_length(reps$indices, 2L)
  expect_equal(dim(reps$var_explained), c(2L, 1L, 2L))
})

test_that(".rajive_bootstrap returns aligned individual score replicates", {
  fx <- make_small_rajive_fixture()
  n_indiv <- ncol(fx$fit$block_decomps[[1L]]$u)

  set.seed(44)
  reps <- rajiveplus:::.rajive_bootstrap(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    B = 2L,
    sample_frac = 0.75,
    keep = "scores",
    score_type = "individual",
    score_block = 1L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  expect_equal(dim(reps$scores), c(12L, n_indiv, 2L))
  expect_true(any(is.finite(reps$scores)))
})

test_that(".rajive_bootstrap is deterministic under set.seed", {
  fx <- make_small_rajive_fixture()

  set.seed(5)
  a <- rajiveplus:::.rajive_bootstrap(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    B = 2L,
    keep = c("joint_rank", "indices"),
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )
  set.seed(5)
  b <- rajiveplus:::.rajive_bootstrap(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    B = 2L,
    keep = c("joint_rank", "indices"),
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  expect_equal(a, b)
})

test_that("assess_stability joint_rank works without fitted reference components", {
  fx <- make_small_rajive_fixture()

  set.seed(51)
  out <- assess_stability(
    blocks = fx$blocks,
    initial_signal_ranks = fx$initial_signal_ranks,
    target = "joint_rank",
    B = 2L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  expect_length(out$rank_distribution, 2L)
  expect_true(inherits(out$rank_table, "table"))
  expect_true(is.na(out$observed_rank))
})

test_that("assess_stability can attach bootstrap replicates", {
  fx <- make_small_rajive_fixture()

  set.seed(6)
  out <- assess_stability(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    target = "components",
    B = 2L,
    return_replicates = TRUE,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  reps <- attr(out, "replicates")
  expect_true(is.data.frame(out))
  expect_true(is.list(reps))
  expect_equal(dim(reps$scores), c(12L, 1L, 2L))
})
