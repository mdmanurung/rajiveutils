library(rajiveplus)

test_that("fortify variance returns tidy long output", {
  fx <- make_extractor_fixture()

  out <- ggplot2::fortify(fx$fit, what = "variance", blocks = fx$blocks)

  expect_true(is.data.frame(out))
  expect_true(all(c("block", "component", "proportion") %in% names(out)))
  expect_equal(nrow(out), length(fx$blocks) * 3L)
  expect_setequal(out$block, seq_along(fx$blocks))
  expect_setequal(out$component, c("Joint", "Indiv", "Resid"))
  expect_true(all(is.finite(out$proportion)))
})

test_that("wide variance extraction keeps legacy component columns", {
  fx <- make_extractor_fixture()

  out <- extract_components(fx$fit, what = "variance", blocks = fx$blocks)

  expect_true(is.data.frame(out))
  expect_named(out, c("Joint", "Indiv", "Resid"))
  expect_equal(nrow(out), length(fx$blocks))
})
