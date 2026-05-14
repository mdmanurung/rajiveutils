make_missing_diagnostics_fixture <- function() {
  set.seed(8601)
  blocks <- list(
    block1 = matrix(rnorm(72), nrow = 12),
    block2 = matrix(rnorm(60), nrow = 12)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[2, 3] <- FALSE
  mask$block1[5, ] <- FALSE
  mask$block2[8, 2] <- FALSE
  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L)
  list(fit = fit, mask = mask)
}

test_that("missingness info exposes summaries and pattern counts", {
  fx <- make_missing_diagnostics_fixture()
  info <- get_missingness_info(fx$fit)

  expect_true(all(c("block_summary", "sample_summary", "counts") %in% names(info)))
  expect_gte(info$counts[["missing_cells"]], 1L)
  expect_gte(info$counts[["whole_sample_block_rows"]], 1L)
})

test_that("missing diagnostics expose variance, support, and convergence", {
  fx <- make_missing_diagnostics_fixture()
  diag <- get_missing_diagnostics(fx$fit)

  expect_true(is.data.frame(diag$variance_explained))
  expect_true(is.data.frame(diag$joint_support))
  expect_true(is.data.frame(diag$orthogonality))
  expect_true(all(is.finite(diag$orthogonality$max_abs_crossprod)))
  expect_true(is.list(diag$convergence))
  expect_true(all(c("Joint", "Indiv", "Resid") %in% names(diag$variance_explained)))
})

test_that("plot_missingness returns a ggplot object", {
  fx <- make_missing_diagnostics_fixture()
  p <- plot_missingness(fx$fit)

  expect_s3_class(p, "ggplot")
})
