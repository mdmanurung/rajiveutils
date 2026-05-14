test_that("native missing mask validation rejects incompatible masks", {
  blocks <- list(
    block1 = matrix(1:12, nrow = 4),
    block2 = matrix(1:8, nrow = 4)
  )
  mask <- list(
    block1 = matrix(TRUE, nrow = 4, ncol = 3),
    block2 = matrix(TRUE, nrow = 3, ncol = 2)
  )

  expect_error(
    rajiveplus:::.normalize_missing_mask(blocks, mask),
    class = "rajiveplus_invalid_input"
  )
})

test_that("native missing validation rejects non-finite observed entries", {
  blocks <- list(
    block1 = matrix(1:12, nrow = 4),
    block2 = matrix(1:8, nrow = 4)
  )
  blocks$block1[2, 1] <- NA_real_
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))

  expect_error(
    rajiveplus:::.normalize_missing_mask(blocks, mask),
    class = "rajiveplus_invalid_input"
  )

  mask$block1[2, 1] <- FALSE
  normalized <- rajiveplus:::.normalize_missing_mask(blocks, mask)
  expect_false(normalized$mask$block1[2, 1])
})

test_that("native missing validation rejects incompatible sample counts", {
  blocks <- list(
    block1 = matrix(1:12, nrow = 4),
    block2 = matrix(1:10, nrow = 5)
  )

  expect_error(
    rajiveplus:::.normalize_missing_mask(blocks),
    class = "rajiveplus_invalid_input"
  )
})

test_that("native missing validation rejects unsupported all-missing structures", {
  blocks <- list(
    block1 = matrix(1:12, nrow = 4),
    block2 = matrix(1:8, nrow = 4)
  )
  rownames(blocks$block1) <- rownames(blocks$block2) <- paste0("s", seq_len(4))
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))

  mask_all_block <- mask
  mask_all_block$block2[,] <- FALSE
  expect_error(
    rajiveplus:::.validate_native_missing_inputs(blocks, mask_all_block),
    class = "rajiveplus_all_missing_block"
  )

  mask_all_sample <- mask
  mask_all_sample$block1[3, ] <- FALSE
  mask_all_sample$block2[3, ] <- FALSE
  expect_error(
    rajiveplus:::.validate_native_missing_inputs(blocks, mask_all_sample),
    regexp = "s3",
    class = "rajiveplus_sample_all_missing"
  )

  mask_all_feature <- mask
  mask_all_feature$block1[, 2] <- FALSE
  expect_error(
    rajiveplus:::.validate_native_missing_inputs(blocks, mask_all_feature),
    class = "rajiveplus_feature_all_missing"
  )
})

test_that("native missing validation rejects rank count mismatches", {
  blocks <- list(
    block1 = matrix(1:12, nrow = 4),
    block2 = matrix(1:8, nrow = 4)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))

  expect_error(
    Rajive(blocks, c(2L), missing = "native", mask = mask, joint_rank = 1L,
           identifiability_norm = "l2"),
    class = "rajiveplus_invalid_input"
  )
  expect_error(
    diagnose_missing_ranks(blocks, candidates = 0:1,
                           initial_signal_ranks = c(2L), mask = mask),
    class = "rajiveplus_invalid_input"
  )
})

test_that("native missing classifier separates common missingness patterns", {
  blocks <- list(
    block1 = matrix(1:20, nrow = 5),
    block2 = matrix(1:15, nrow = 5)
  )
  mask <- lapply(blocks, function(x) matrix(TRUE, nrow(x), ncol(x)))
  mask$block1[2, ] <- FALSE
  mask$block1[4, 3] <- FALSE
  mask$block2[, 2] <- FALSE

  patterns <- rajiveplus:::.classify_missingness(blocks, mask)

  expect_true(patterns$sample_block_rows$whole_row_missing[patterns$sample_block_rows$block == "block1" &
                                                            patterns$sample_block_rows$sample == 2])
  expect_true(patterns$cells$scattered_missing[patterns$cells$block == "block1" &
                                                 patterns$cells$row == 4 &
                                                 patterns$cells$col == 3])
  expect_true(patterns$features$all_missing[patterns$features$block == "block2" &
                                              patterns$features$feature == 2])
  expect_equal(patterns$counts[["whole_sample_block_rows"]], 1L)
  expect_equal(patterns$counts[["all_missing_features"]], 1L)
})
