test_that("rajive_ci returns percentile loading intervals", {
  fx <- make_small_rajive_fixture()

  set.seed(11)
  out <- rajive_ci(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    target = "loadings",
    method = "percentile",
    B = 3L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  expected <- c("target", "block", "component", "feature", "sample",
                "estimate", "lower", "upper", "level", "method",
                "n_replicates")
  expect_named(out, expected)
  expect_true(all(out$target == "loadings"))
  expect_true(all(out$method == "percentile"))
  expect_equal(nrow(out), sum(vapply(fx$blocks, ncol, integer(1))))
})

test_that("rajive_ci returns basic variance-explained intervals", {
  fx <- make_small_rajive_fixture()

  set.seed(12)
  out <- rajive_ci(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    target = "var_explained",
    method = "basic",
    B = 3L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  expect_true(all(out$target == "var_explained"))
  expect_true(all(out$method == "basic"))
  expect_equal(nrow(out), length(fx$blocks))
  expect_true(all(is.finite(out$estimate)))
})

test_that("rajive_ci computes BCa intervals using jackknife refits", {
  fx <- make_small_rajive_fixture()

  set.seed(13)
  out <- rajive_ci(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    target = "joint_rank",
    method = "bca",
    B = 3L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )

  expect_equal(out$target, "joint_rank")
  expect_equal(out$method, "bca")
  expect_equal(out$n_replicates, 3L)
  expect_true(is.finite(out$lower))
  expect_true(is.finite(out$upper))
})

test_that("jackknife refits inherit fitted identifiability_norm", {
  fx <- make_small_rajive_fixture()
  fx$fit$joint_rank_sel$identifiability_norm <- "l1"

  seen <- character()
  fake_fit <- function(b_list, initial_signal_ranks,
                       identifiability_norm = NULL, ...) {
    seen <<- c(seen, identifiability_norm)
    structure(
      list(
        joint_rank = 1L,
        joint_scores = matrix(seq_len(nrow(b_list[[1L]])), ncol = 1L)
      ),
      class = "rajive"
    )
  }

  testthat::with_mocked_bindings(
    rajiveplus:::.rajive_jackknife(
      fx$fit,
      fx$blocks,
      fx$initial_signal_ranks,
      target = "joint_rank",
      joint_rank = 1L,
      n_wedin_samples = NA,
      n_rand_dir_samples = NA
    ),
    Rajive = fake_fit,
    .package = "rajiveplus"
  )

  expect_equal(seen, rep("l1", nrow(fx$blocks[[1L]])))
})

test_that("rajive_ci rejects BCa for sample-specific scores", {
  fx <- make_small_rajive_fixture()

  expect_error(
    rajive_ci(
      fx$fit, fx$blocks, fx$initial_signal_ranks,
      target = "scores",
      method = "bca",
      B = 3L
    ),
    regexp = "BCa"
  )
})

test_that("extract_components can merge loading confidence intervals", {
  fx <- make_small_rajive_fixture()

  set.seed(14)
  ci <- rajive_ci(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    target = "loadings",
    method = "percentile",
    B = 3L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )
  out <- extract_components(fx$fit, what = "loadings", format = "long", ci = ci)

  expect_true(all(c("lower", "upper", "level", "ci_method") %in% names(out)))
  expect_equal(nrow(out), nrow(ci))
})

test_that("extract_components merges loading confidence intervals by feature key", {
  fx <- make_small_rajive_fixture()

  set.seed(15)
  ci <- rajive_ci(
    fx$fit, fx$blocks, fx$initial_signal_ranks,
    target = "loadings",
    method = "percentile",
    B = 3L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )
  ci <- ci[rev(seq_len(nrow(ci))), , drop = FALSE]
  out <- extract_components(fx$fit, what = "loadings", format = "long", ci = ci)

  key <- paste(out$block, out$component, out$feature)
  ci_key <- paste(
    sub("^block", "", ci$block),
    ci$component,
    sub("^feature", "", ci$feature)
  )
  expect_equal(out$lower, ci$lower[match(key, ci_key)])
  expect_equal(out$upper, ci$upper[match(key, ci_key)])
})

test_that("extract_components merges shuffled loading intervals with named blocks", {
  set.seed(16)
  Y <- ajive.data.sim(K = 2, rankJ = 1, rankA = c(1, 1),
                      n = 12, pks = c(8, 6), dist.type = 1)
  blocks <- Y$sim_data
  names(blocks) <- c("rna", "protein")
  colnames(blocks$rna) <- paste0("gene", seq_len(ncol(blocks$rna)))
  colnames(blocks$protein) <- paste0("protein", seq_len(ncol(blocks$protein)))
  fit <- Rajive(blocks, c(2L, 2L),
                joint_rank = 1L,
                n_wedin_samples = NA,
                n_rand_dir_samples = NA,
                num_cores = 1L)

  set.seed(17)
  ci <- rajive_ci(
    fit, blocks, c(2L, 2L),
    target = "loadings",
    method = "percentile",
    B = 3L,
    n_wedin_samples = NA,
    n_rand_dir_samples = NA,
    joint_rank = 1L
  )
  ci <- ci[rev(seq_len(nrow(ci))), , drop = FALSE]
  out <- extract_components(fit, what = "loadings", format = "long", ci = ci)

  expect_equal(out$lower, ci$lower[match(out$loading, ci$estimate)])
  expect_equal(out$upper, ci$upper[match(out$loading, ci$estimate)])
})
