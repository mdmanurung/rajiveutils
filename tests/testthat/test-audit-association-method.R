library(rajiveplus)

make_audit_assoc_fit <- function() {
  n <- 8L
  joint_scores <- cbind(
    comp1 = seq_len(n),
    comp2 = c(2, 1, 4, 3, 6, 5, 8, 7)
  )
  indiv <- list(
    u = matrix(seq_len(n), ncol = 1L),
    d = 1,
    v = matrix(1, nrow = 1L, ncol = 1L)
  )
  joint <- list(
    u = joint_scores,
    d = c(2, 1),
    v = diag(2)
  )
  noise <- matrix(0, nrow = n, ncol = 1L)

  structure(
    list(
      joint_scores = joint_scores,
      joint_rank = 2L,
      block_decomps = list(indiv, joint, noise, indiv, joint, noise)
    ),
    class = "rajive"
  )
}

test_that("categorical association rejects unsupported method labels", {
  fit <- make_audit_assoc_fit()
  metadata <- data.frame(group = rep(c("a", "b"), each = 4L))

  expect_error(
    associate_components(
      fit,
      metadata,
      variable = "group",
      mode = "categorical",
      method = "anova"
    ),
    class = "rajiveplus_invalid_input"
  )
})

test_that("categorical and batch associations report the executed kruskal test", {
  fit <- make_audit_assoc_fit()
  metadata <- data.frame(group = rep(c("a", "b"), each = 4L))

  cat_default <- suppressMessages(
    associate_components(fit, metadata, variable = "group", mode = "categorical")
  )
  cat_explicit <- suppressMessages(
    associate_components(fit, metadata, variable = "group",
                         mode = "categorical", method = "kruskal")
  )
  batch_explicit <- suppressMessages(
    associate_components(fit, metadata, variable = "group",
                         mode = "batch", method = "kruskal")
  )

  expect_true(all(cat_default$method == "kruskal"))
  expect_true(all(cat_explicit$method == "kruskal"))
  expect_true(all(batch_explicit$method == "kruskal"))
})

test_that("categorical wrapper and all-component wrapper preserve kruskal labels", {
  fit <- make_audit_assoc_fit()
  metadata <- data.frame(group = rep(c("a", "b"), each = 4L))

  wrapped <- suppressMessages(
    associate_scores_categorical(fit, metadata, variable = "group")
  )
  all_components <- suppressMessages(
    associate_all_components(
      fit,
      metadata,
      variable = "group",
      mode = "categorical",
      include = "both",
      blocks_to_include = 1L
    )
  )

  expect_true(all(wrapped$method == "kruskal"))
  expect_true(all(all_components$method == "kruskal"))
})
