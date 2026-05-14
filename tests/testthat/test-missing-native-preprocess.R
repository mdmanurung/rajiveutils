test_that("observed-entry preprocessing centers only observed values", {
  X <- matrix(c(1, 2, 1000,
                4, 5, 6,
                7, 8, 9), nrow = 3, byrow = TRUE)
  blocks <- list(block1 = X)
  mask <- list(block1 = matrix(TRUE, nrow = 3, ncol = 3))
  mask$block1[1, 3] <- FALSE

  prep <- rajiveplus:::.center_scale_observed(blocks, mask,
                                              center = TRUE,
                                              scale = FALSE,
                                              normalize = FALSE)

  observed_means <- colSums(prep$blocks$block1 * prep$mask$block1, na.rm = TRUE) /
    colSums(prep$mask$block1)
  expect_equal(observed_means, c(0, 0, 0), tolerance = 1e-12)
  expect_true(is.na(prep$blocks$block1[1, 3]))

  X2 <- X
  X2[1, 3] <- -999999
  prep2 <- rajiveplus:::.center_scale_observed(list(block1 = X2), mask,
                                               center = TRUE,
                                               scale = FALSE,
                                               normalize = FALSE)
  expect_equal(prep$center$block1, prep2$center$block1)
  expect_equal(prep$blocks$block1[mask$block1], prep2$blocks$block1[mask$block1])
})

test_that("observed-entry preprocessing can normalize observed Frobenius norm", {
  X <- matrix(c(1, 2, NA,
                4, 5, 6,
                7, 8, 9), nrow = 3, byrow = TRUE)
  mask <- list(block1 = is.finite(X))
  prep <- rajiveplus:::.center_scale_observed(list(block1 = X), mask,
                                              center = TRUE,
                                              scale = FALSE,
                                              normalize = TRUE)

  observed_norm <- sqrt(sum(prep$blocks$block1[mask$block1]^2))
  expect_equal(observed_norm, 1, tolerance = 1e-12)
  expect_true(is.finite(prep$frob_norm$block1))
})

test_that("backtransform restores observed entries after centering and scaling", {
  X <- matrix(c(1, 2, NA,
                4, 5, 6,
                7, 8, 9), nrow = 3, byrow = TRUE)
  mask <- list(block1 = is.finite(X))
  prep <- rajiveplus:::.center_scale_observed(list(block1 = X), mask,
                                              center = TRUE,
                                              scale = TRUE,
                                              normalize = TRUE)

  restored <- rajiveplus:::.backtransform_reconstruction(prep$blocks, prep)
  expect_equal(restored$block1[mask$block1], X[mask$block1], tolerance = 1e-10)
  expect_true(is.na(restored$block1[!mask$block1]))
})

test_that("all-finite preprocessing records complete observed counts", {
  X <- matrix(rnorm(20), nrow = 5)
  mask <- list(block1 = matrix(TRUE, nrow = 5, ncol = 4))
  prep <- rajiveplus:::.center_scale_observed(list(block1 = X), mask,
                                              center = TRUE,
                                              scale = FALSE,
                                              normalize = TRUE)

  expect_equal(prep$observed_counts$block1, rep(5L, 4L))
  expect_equal(sqrt(sum(prep$blocks$block1^2)), 1, tolerance = 1e-12)
})
