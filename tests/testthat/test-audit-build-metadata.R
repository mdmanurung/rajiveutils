test_that("pkgdown article index matches package vignette surface", {
  root <- getwd()
  if (!file.exists(file.path(root, "_pkgdown.yml"))) {
    candidates <- c(file.path(root, "..", ".."), file.path(root, ".."))
    candidates <- candidates[file.exists(file.path(candidates, "_pkgdown.yml"))]
    skip_if(
      length(candidates) == 0L,
      "pkgdown metadata is not available in this check environment"
    )
    root <- candidates[[1L]]
  }

  rmd_paths <- list.files(file.path(root, "vignettes"),
                          pattern = "\\.Rmd$", full.names = FALSE)
  vignette_ids <- tools::file_path_sans_ext(basename(rmd_paths))

  expect_false("benchmarking_heavy" %in% vignette_ids)

  indexed <- trimws(readLines(file.path(root, "_pkgdown.yml"), warn = FALSE))
  indexed <- indexed[grepl("^- [A-Za-z0-9_]+$", indexed)]
  indexed <- sub("^- ", "", indexed)

  missing <- setdiff(vignette_ids, indexed)
  expect_equal(missing, character(0))
})
