test_that("associate_components default output keeps legacy columns", {
  set.seed(21)
  scores <- matrix(rnorm(30), ncol = 1)
  fit <- structure(list(joint_scores = scores, joint_rank = 1L), class = "rajive")
  metadata <- data.frame(marker = rnorm(30))

  out <- suppressMessages(
    associate_components(fit, metadata, variable = "marker",
                         mode = "continuous")
  )

  expect_named(out, c("variable", "component", "stat", "p_value",
                      "p_adj", "method"))
})

test_that("associate_components bootstrap propagation adds uncertainty columns", {
  set.seed(22)
  n <- 30L
  B <- 6L
  score <- rnorm(n)
  metadata <- data.frame(marker = score + rnorm(n, sd = 0.05))
  fit <- structure(list(joint_scores = matrix(score, ncol = 1),
                        joint_rank = 1L),
                   class = "rajive")
  reps <- list(scores = array(rep(score, B), dim = c(n, 1L, B)))
  for (b in seq_len(B)) {
    reps$scores[, 1L, b] <- score + rnorm(n, sd = 0.02)
  }

  out <- suppressMessages(
    associate_components(
      fit, metadata, variable = "marker",
      mode = "continuous",
      propagate_uncertainty = "bootstrap",
      alpha_stability = 0.05,
      replicates = reps
    )
  )

  expect_true(all(c("stability", "effect_lo", "effect_hi", "p_median",
                    "p_adj_median", "n_boot_valid") %in% names(out)))
  expect_gte(out$stability, 0.9)
  expect_gt(out$effect_lo, 0)
  expect_lt(out$p_median, 0.05)
  expect_equal(out$n_boot_valid, B)
})

test_that("associate_components bootstrap propagation can build replicates", {
  fx <- make_small_rajive_fixture(seed = 23L)
  score <- fx$fit$joint_scores[, 1L]
  metadata <- data.frame(marker = score + rnorm(length(score), sd = 0.05))

  set.seed(24)
  out <- suppressMessages(
    associate_components(
      fx$fit, metadata, variable = "marker",
      mode = "continuous",
      propagate_uncertainty = "bootstrap",
      blocks = fx$blocks,
      initial_signal_ranks = fx$initial_signal_ranks,
      B = 2L,
      n_wedin_samples = NA,
      n_rand_dir_samples = NA,
      joint_rank = 1L
    )
  )

  expect_true(all(c("stability", "effect_lo", "effect_hi", "p_median",
                    "p_adj_median", "n_boot_valid") %in% names(out)))
  expect_equal(nrow(out), 1L)
  expect_gte(out$n_boot_valid, 0L)
})

test_that("associate_components bootstrap propagation supports individual scores", {
  fx <- make_small_rajive_fixture(seed = 25L)
  score <- fx$fit$block_decomps[[1L]]$u[, 1L]
  metadata <- data.frame(marker = score + rnorm(length(score), sd = 0.05))

  set.seed(26)
  out <- suppressMessages(
    associate_components(
      fx$fit, metadata, variable = "marker",
      mode = "continuous",
      type = "individual",
      block = 1L,
      propagate_uncertainty = "bootstrap",
      blocks = fx$blocks,
      initial_signal_ranks = fx$initial_signal_ranks,
      B = 2L,
      n_wedin_samples = NA,
      n_rand_dir_samples = NA,
      joint_rank = 1L
    )
  )

  expect_true(all(c("stability", "effect_lo", "effect_hi", "p_median",
                    "p_adj_median", "n_boot_valid") %in% names(out)))
  expect_equal(nrow(out), ncol(fx$fit$block_decomps[[1L]]$u))
  expect_true(any(out$n_boot_valid > 0L))
})

test_that("associate_all_components returns requested row identities", {
  fx <- make_small_rajive_fixture(seed = 27L)
  metadata <- data.frame(marker = rnorm(nrow(fx$fit$joint_scores)))

  joint <- suppressMessages(
    associate_all_components(fx$fit, metadata, variable = "marker",
                             mode = "continuous", include = "joint")
  )
  indiv <- suppressMessages(
    associate_all_components(fx$fit, metadata, variable = "marker",
                             mode = "continuous", include = "individual",
                             blocks_to_include = 1L)
  )
  both <- suppressMessages(
    associate_all_components(fx$fit, metadata, variable = "marker",
                             mode = "continuous", include = "both",
                             blocks_to_include = 1L)
  )

  expect_true(all(joint$source == "joint"))
  expect_true(all(is.na(joint$block)))
  expect_true(all(indiv$source == "individual"))
  expect_equal(unique(indiv$block), 1L)
  expect_equal(nrow(both), nrow(joint) + nrow(indiv))
  expect_named(
    both,
    c("source", "block", "variable", "component", "stat", "p_value",
      "p_adj", "method")
  )
})

test_that("associate_all_components bootstraps joint and individual targets together", {
  fx <- make_small_rajive_fixture(seed = 28L)
  metadata <- data.frame(marker = fx$fit$joint_scores[, 1L] + rnorm(12, sd = 0.1))

  set.seed(29)
  out <- suppressMessages(
    associate_all_components(
      fx$fit, metadata, variable = "marker",
      mode = "continuous",
      include = "both",
      blocks_to_include = 1L,
      propagate_uncertainty = "bootstrap",
      blocks = fx$blocks,
      initial_signal_ranks = fx$initial_signal_ranks,
      B = 2L,
      n_wedin_samples = NA,
      n_rand_dir_samples = NA,
      joint_rank = 1L
    )
  )

  expect_true(all(c("stability", "effect_lo", "effect_hi", "p_median",
                    "p_adj_median", "n_boot_valid") %in% names(out)))
  expect_true(all(c("joint", "individual") %in% out$source))
  expect_equal(out$p_adj, stats::p.adjust(out$p_value, method = "BH"))
  expect_equal(out$p_adj_median, stats::p.adjust(out$p_median, method = "BH"))
  expect_true(any(out$n_boot_valid > 0L))
})
