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

test_that(".rajive_bootstrap rank-only payload matches full refit ranks and indices", {
  fx <- make_small_rajive_fixture()
  cluster <- rep(letters[1:4], each = 3L)

  set.seed(51)
  rank_only <- rajiveplus:::.rajive_bootstrap(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    B = 3L,
    keep = c("joint_rank", "indices"),
    cluster = cluster,
    resample = "cluster",
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  set.seed(51)
  full <- rajiveplus:::.rajive_bootstrap(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    B = 3L,
    keep = c("joint_rank", "indices", "scores"),
    cluster = cluster,
    resample = "cluster",
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  expect_equal(rank_only$joint_rank, full$joint_rank)
  expect_equal(rank_only$indices, full$indices)
})

test_that("bootstrap refits use num_cores for outer replicates and serial inner Rajive calls", {
  n <- 8L
  blocks <- list(
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L)),
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L))
  )

  seen_single <- integer()
  fake_rank <- function(b_list, initial_signal_ranks, num_cores = 1L, ...) {
    seen_single <<- c(seen_single, as.integer(num_cores))
    structure(list(joint_rank = 1L), class = "rajive_rank_only")
  }
  fake_parallel <- function(X, FUN, num_cores = 1L, ...) {
    lapply(X, FUN, ...)
  }

  set.seed(7)
  testthat::with_mocked_bindings(
    rajiveplus:::.rajive_bootstrap(
      ajive_output = NULL,
      blocks = blocks,
      initial_signal_ranks = c(1L, 1L),
      B = 2L,
      keep = "joint_rank",
      num_cores = 5L
    ),
    .rajive_parallel_lapply = fake_parallel,
    .Rajive_rank_only = fake_rank,
    .package = "rajiveplus"
  )
  expect_equal(seen_single, c(1L, 1L))

  ref_scores <- matrix(seq_len(n), ncol = 1L)
  ajive_output <- structure(list(joint_scores = ref_scores), class = "rajive")
  targets <- data.frame(
    source = "joint",
    block = NA_integer_,
    key = "joint",
    stringsAsFactors = FALSE
  )

  seen_multi <- integer()
  fake_scores <- function(b_list, initial_signal_ranks, num_cores = 1L, ...) {
    seen_multi <<- c(seen_multi, as.integer(num_cores))
    list(joint_scores = matrix(b_list[[1L]][, 1L], ncol = 1L), joint_rank = 1L)
  }

  set.seed(8)
  testthat::with_mocked_bindings(
    rajiveplus:::.rajive_bootstrap_scores_multi(
      ajive_output = ajive_output,
      blocks = blocks,
      initial_signal_ranks = c(1L, 1L),
      targets = targets,
      B = 2L,
      num_cores = 6L
    ),
    .rajive_parallel_lapply = fake_parallel,
    Rajive = fake_scores,
    .package = "rajiveplus"
  )
  expect_equal(seen_multi, c(1L, 1L))
})

test_that("bootstrap refits are parallelized at replicate level when num_cores > 1", {
  n <- 8L
  blocks <- list(
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L)),
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L))
  )

  seen_parallel <- FALSE
  seen_inner_cores <- integer()
  fake_parallel <- function(X, FUN, num_cores = 1L, ...) {
    seen_parallel <<- TRUE
    expect_equal(num_cores, 2L)
    lapply(X, FUN, ...)
  }
  fake_rank_only <- function(b_list, initial_signal_ranks, num_cores = 1L, ...) {
    seen_inner_cores <<- c(seen_inner_cores, as.integer(num_cores))
    structure(list(joint_rank = 1L), class = "rajive_rank_only")
  }

  set.seed(71)
  testthat::with_mocked_bindings(
    rajiveplus:::.rajive_bootstrap(
      ajive_output = NULL,
      blocks = blocks,
      initial_signal_ranks = c(1L, 1L),
      B = 3L,
      keep = "joint_rank",
      num_cores = 2L
    ),
    .rajive_parallel_lapply = fake_parallel,
    .Rajive_rank_only = fake_rank_only,
    .package = "rajiveplus"
  )

  expect_true(seen_parallel)
  expect_equal(seen_inner_cores, rep(1L, 3L))
})

test_that("bootstrap refits inherit fitted identifiability_norm", {
  n <- 8L
  blocks <- list(
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L)),
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L))
  )
  ref_scores <- matrix(seq_len(n), ncol = 1L)
  ajive_output <- structure(
    list(
      joint_scores = ref_scores,
      joint_rank_sel = list(identifiability_norm = "l1")
    ),
    class = "rajive"
  )

  seen <- character()
  fake_full <- function(b_list, initial_signal_ranks,
                        identifiability_norm = NULL, ...) {
    seen <<- c(seen, identifiability_norm)
    list(joint_scores = matrix(b_list[[1L]][, 1L], ncol = 1L), joint_rank = 1L)
  }

  set.seed(72)
  testthat::with_mocked_bindings(
    rajiveplus:::.rajive_bootstrap(
      ajive_output = ajive_output,
      blocks = blocks,
      initial_signal_ranks = c(1L, 1L),
      B = 2L,
      keep = "scores"
    ),
    Rajive = fake_full,
    .package = "rajiveplus"
  )
  expect_equal(seen, c("l1", "l1"))

  seen <- character()
  set.seed(73)
  testthat::with_mocked_bindings(
    rajiveplus:::.rajive_bootstrap(
      ajive_output = ajive_output,
      blocks = blocks,
      initial_signal_ranks = c(1L, 1L),
      B = 2L,
      keep = "scores",
      identifiability_norm = "l2"
    ),
    Rajive = fake_full,
    .package = "rajiveplus"
  )
  expect_equal(seen, c("l2", "l2"))
})

test_that("rank-only bootstrap refits default identifiability_norm to L2 quietly", {
  n <- 8L
  blocks <- list(
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L)),
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L))
  )

  seen <- character()
  fake_rank_only <- function(b_list, initial_signal_ranks,
                             identifiability_norm = NULL, ...) {
    seen <<- c(seen, identifiability_norm)
    structure(list(joint_rank = 1L), class = "rajive_rank_only")
  }

  set.seed(74)
  expect_no_message(
    testthat::with_mocked_bindings(
      rajiveplus:::.rajive_bootstrap(
        ajive_output = NULL,
        blocks = blocks,
        initial_signal_ranks = c(1L, 1L),
        B = 2L,
        keep = "joint_rank"
      ),
      .Rajive_rank_only = fake_rank_only,
      .package = "rajiveplus"
    )
  )
  expect_equal(seen, c("l2", "l2"))
})

test_that("multi-score bootstrap refits inherit fitted identifiability_norm", {
  n <- 8L
  blocks <- list(
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L)),
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L))
  )
  ajive_output <- structure(
    list(
      joint_scores = matrix(seq_len(n), ncol = 1L),
      joint_rank_sel = list(identifiability_norm = "l1")
    ),
    class = "rajive"
  )
  targets <- data.frame(
    source = "joint",
    block = NA_integer_,
    key = "joint",
    stringsAsFactors = FALSE
  )

  seen <- character()
  fake_full <- function(b_list, initial_signal_ranks,
                        identifiability_norm = NULL, ...) {
    seen <<- c(seen, identifiability_norm)
    list(joint_scores = matrix(b_list[[1L]][, 1L], ncol = 1L), joint_rank = 1L)
  }

  set.seed(75)
  testthat::with_mocked_bindings(
    rajiveplus:::.rajive_bootstrap_scores_multi(
      ajive_output = ajive_output,
      blocks = blocks,
      initial_signal_ranks = c(1L, 1L),
      targets = targets,
      B = 2L
    ),
    Rajive = fake_full,
    .package = "rajiveplus"
  )
  expect_equal(seen, c("l1", "l1"))
})

test_that(".rajive_bootstrap uses rank-only refits only for rank-only payloads", {
  n <- 8L
  blocks <- list(
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L)),
    cbind(seq_len(n), matrix(rnorm(n * 2L), n, 2L))
  )
  ref_scores <- matrix(seq_len(n), ncol = 1L)
  ajive_output <- structure(list(joint_scores = ref_scores), class = "rajive")

  full_calls <- 0L
  rank_only_calls <- 0L
  fake_full <- function(b_list, initial_signal_ranks, ...) {
    full_calls <<- full_calls + 1L
    list(joint_scores = matrix(b_list[[1L]][, 1L], ncol = 1L), joint_rank = 1L)
  }
  fake_rank_only <- function(b_list, initial_signal_ranks, ...) {
    rank_only_calls <<- rank_only_calls + 1L
    structure(list(joint_scores = matrix(b_list[[1L]][, 1L], ncol = 1L),
                   joint_rank = 1L),
              class = "rajive_rank_only")
  }

  set.seed(9)
  testthat::with_mocked_bindings(
    rajiveplus:::.rajive_bootstrap(
      ajive_output = NULL,
      blocks = blocks,
      initial_signal_ranks = c(1L, 1L),
      B = 2L,
      keep = c("joint_rank", "indices")
    ),
    Rajive = fake_full,
    .Rajive_rank_only = fake_rank_only,
    .package = "rajiveplus"
  )
  expect_equal(rank_only_calls, 2L)
  expect_equal(full_calls, 0L)

  set.seed(10)
  testthat::with_mocked_bindings(
    rajiveplus:::.rajive_bootstrap(
      ajive_output = ajive_output,
      blocks = blocks,
      initial_signal_ranks = c(1L, 1L),
      B = 2L,
      keep = "scores"
    ),
    Rajive = fake_full,
    .Rajive_rank_only = fake_rank_only,
    .package = "rajiveplus"
  )
  expect_equal(rank_only_calls, 2L)
  expect_equal(full_calls, 2L)
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
