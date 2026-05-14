make_missing_estimability_fixture <- function() {
  set.seed(8301)
  blocks <- list(
    block1 = matrix(rnorm(60), nrow = 10),
    block2 = matrix(rnorm(50), nrow = 10)
  )
  rownames(blocks$block1) <- rownames(blocks$block2) <- paste0("s", seq_len(10))
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[1, 1] <- FALSE
  mask$block2[4, ] <- FALSE
  fit <- rajiveplus:::.Rajive_incomplete(blocks, c(2L, 2L),
                                         joint_rank = 1L, mask = mask)
  list(fit = fit, mask = mask)
}

make_weak_support_fixture <- function() {
  set.seed(8302)
  blocks <- list(block1 = matrix(rnorm(60), nrow = 10),
                 block2 = matrix(rnorm(50), nrow = 10))
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[3:10, 2] <- FALSE
  fit <- Rajive(blocks, c(2L, 2L), missing = "native", mask = mask,
                joint_rank = 1L,
                missing_control = rajive_missing_control(min_component_support = 3L))
  list(fit = fit, mask = mask)
}

test_that("estimability labels observed and scattered missing cells", {
  fx <- make_missing_estimability_fixture()
  est <- get_estimability(fx$fit, block = "block1")

  expect_equal(est$label[est$row == 2 & est$col == 1 &
                           est$component == "joint_individual"],
               "observed")
  expect_equal(est$label[est$row == 1 & est$col == 1 &
                           est$component == "joint_individual"],
               "within_block_estimable")
})

test_that("estimability labels sample-block rows honestly", {
  fx <- make_missing_estimability_fixture()
  est <- get_estimability(fx$fit, block = "block2")

  expect_equal(unique(est$label[est$row == 4 & est$component == "joint"]),
               "cross_block_joint_estimable")
  expect_equal(unique(est$label[est$row == 4 & est$component == "individual"]),
               "not_identifiable")
})

test_that("weak estimability reflects min_component_support and warns with a class", {
  fx <- make_weak_support_fixture()
  est <- get_estimability(fx$fit, block = "block1")
  expect_true(any(est$label == "weakly_estimable"))
  expect_warning(get_missing_diagnostics(fx$fit, warn_weak = TRUE),
                 class = "rajiveplus_weak_estimability")
})
