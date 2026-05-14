test_that("native missing union-sample vignette renders", {
  skip_if_not_installed("rmarkdown")

  vignette_path <- test_path("../../vignettes/native_missing_union.Rmd")
  if (!file.exists(vignette_path)) {
    vignette_path <- file.path(getwd(), "vignettes", "native_missing_union.Rmd")
  }
  skip_if_not(file.exists(vignette_path),
              "source vignette is not available in this check environment")

  out_dir <- tempfile("rajiveplus-native-missing-union-")
  dir.create(out_dir)
  rendered <- rmarkdown::render(
    input = vignette_path,
    output_dir = out_dir,
    intermediates_dir = out_dir,
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )

  expect_true(file.exists(rendered))
})
