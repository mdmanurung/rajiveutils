library(rajiveplus)

test_that("assess_stability rejects unsupported permutation method", {
  fx <- make_small_rajive_fixture(seed = 9101L)

  expect_error(
    assess_stability(
      ajive_output = fx$fit,
      blocks = fx$blocks,
      initial_signal_ranks = fx$initial_signal_ranks,
      target = "joint_rank",
      method = "permutation",
      B = 2L,
      n_perm = 7L,
      joint_rank = 1L,
      n_wedin_samples = NA,
      n_rand_dir_samples = NA
    ),
    class = "rajiveplus_invalid_input"
  )
})

test_that("assess_stability bootstrap still uses B for joint-rank summaries", {
  fx <- make_small_rajive_fixture(seed = 9102L)

  set.seed(9103L)
  out <- assess_stability(
    ajive_output = fx$fit,
    blocks = fx$blocks,
    initial_signal_ranks = fx$initial_signal_ranks,
    target = "joint_rank",
    method = "bootstrap",
    B = 3L,
    joint_rank = 1L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA
  )

  expect_equal(length(out$rank_distribution), 3L)
})
