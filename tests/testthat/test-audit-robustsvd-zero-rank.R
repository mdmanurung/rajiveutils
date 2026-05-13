library(rajiveplus)

test_that("RobRSVD.all returns a zero-column decomposition for nrank = 0", {
  x <- matrix(rnorm(20), nrow = 5L, ncol = 4L)

  out <- rajiveplus:::RobRSVD.all(x, nrank = 0L)

  expect_equal(out$d, numeric(0))
  expect_equal(dim(out$u), c(nrow(x), 0L))
  expect_equal(dim(out$v), c(ncol(x), 0L))
})
