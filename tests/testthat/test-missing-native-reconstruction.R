make_missing_reconstruction_fixture <- function() {
  set.seed(8501)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[2, 3] <- FALSE
  mask$block2[4, ] <- FALSE
  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L,
                missing_control = rajive_missing_control(center = TRUE,
                                                         scale = TRUE,
                                                         normalize = TRUE))
  list(fit = fit, blocks = blocks, mask = mask)
}

test_that("reconstruction accessors return finite matrices with original dimensions", {
  fx <- make_missing_reconstruction_fixture()
  recon <- get_reconstructed_blocks(fx$fit, type = "joint")

  expect_equal(length(recon), 2L)
  expect_equal(dim(recon$block1), dim(fx$blocks$block1))
  expect_true(all(is.finite(recon$block1)))
})

test_that("joint_individual reconstruction excludes non-identifiable individual rows", {
  fx <- make_missing_reconstruction_fixture()
  joint <- get_reconstructed_blocks(fx$fit, type = "joint", scale = "standardized")
  both <- get_reconstructed_blocks(fx$fit, type = "joint_individual",
                                   scale = "standardized")

  expect_equal(both$block2[4, ], joint$block2[4, ], tolerance = 1e-10)
  est <- get_estimability(fx$fit, block = "block2")
  expect_equal(unique(est$label[est$row == 4 & est$component == "individual"]),
               "not_identifiable")
})

test_that("observed entries can be returned as observed values", {
  fx <- make_missing_reconstruction_fixture()
  recon <- get_reconstructed_blocks(fx$fit, type = "joint_individual",
                                    observed = "observed")
  expect_equal(recon$block1[fx$mask$block1],
               fx$blocks$block1[fx$mask$block1],
               tolerance = 1e-10)
})

test_that("original-scale reconstruction matches standardized transform", {
  fx <- make_missing_reconstruction_fixture()
  std <- get_reconstructed_blocks(fx$fit, type = "joint_individual",
                                  scale = "standardized")
  orig <- get_reconstructed_blocks(fx$fit, type = "joint_individual",
                                   scale = "original")
  manual <- rajiveplus:::.backtransform_reconstruction(std,
                                                       fx$fit$missing$preprocess)
  expect_equal(orig$block1, manual$block1, tolerance = 1e-10)
})

test_that("reconstruction provenance is stored and returned", {
  fx <- make_missing_reconstruction_fixture()
  recon <- get_reconstructed_blocks(fx$fit, type = "joint")
  provenance <- attr(recon, "reconstruction_provenance", exact = TRUE)

  expect_true(is.data.frame(fx$fit$missing$reconstruction_provenance))
  expect_true(is.data.frame(provenance))
  expect_true(all(provenance$component == "joint"))
  expect_true(all(c("block", "row", "col", "component", "label") %in%
                    names(provenance)))
})
