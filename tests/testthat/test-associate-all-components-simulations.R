benchmark_associate_all_components_scenario <- function(seed, metadata_source,
                                                        B = 6L) {
  set.seed(seed)
  Y <- ajive.data.sim(K = 2, rankJ = 1, rankA = c(1, 1),
                      n = 24, pks = c(16, 14), dist.type = 1)
  fit <- Rajive(Y$sim_data, c(2L, 2L),
                joint_rank = 1L,
                n_wedin_samples = NA,
                n_rand_dir_samples = NA,
                num_cores = 1L)

  marker <- switch(
    metadata_source,
    joint = fit$joint_scores[, 1L] + stats::rnorm(24, sd = 0.15),
    individual = fit$block_decomps[[1L]]$u[, 1L] + stats::rnorm(24, sd = 0.15),
    null = stats::rnorm(24)
  )

  elapsed <- system.time({
    out <- associate_all_components(
      fit,
      data.frame(marker = marker),
      variable = "marker",
      mode = "continuous",
      include = "both",
      propagate_uncertainty = "bootstrap",
      blocks = Y$sim_data,
      initial_signal_ranks = c(2L, 2L),
      B = B,
      n_wedin_samples = NA,
      n_rand_dir_samples = NA,
      joint_rank = 1L
    )
  })

  list(
    result = out,
    metrics = data.frame(
      scenario = metadata_source,
      B = B,
      elapsed_sec = unname(elapsed[["elapsed"]]),
      refit_valid_fraction = mean(out$n_boot_valid / B, na.rm = TRUE),
      min_p_value = min(out$p_value, na.rm = TRUE),
      max_stability = max(out$stability, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  )
}

test_that("associate_all_components slow simulation benchmarks collect diagnostics", {
  skip_if_not_slow()

  scenarios <- lapply(
    c("joint", "individual", "null"),
    function(src) benchmark_associate_all_components_scenario(
      seed = switch(src, joint = 101L, individual = 102L, null = 103L),
      metadata_source = src
    )
  )
  metrics <- do.call(rbind, lapply(scenarios, `[[`, "metrics"))

  expect_equal(metrics$scenario, c("joint", "individual", "null"))
  expect_true(all(is.finite(metrics$elapsed_sec)))
  expect_true(all(metrics$refit_valid_fraction >= 0 &
                    metrics$refit_valid_fraction <= 1))
  expect_true(all(is.finite(metrics$min_p_value)))
  expect_true(all(is.finite(metrics$max_stability)))
})
