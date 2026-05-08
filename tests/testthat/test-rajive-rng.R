# Reproducibility tests for parallel RNG (W-H2).

test_that("get_random_direction_bound_robustH is reproducible under L'Ecuyer", {
  skip_on_cran()

  d <- c(8, 7)
  n <- 30

  a <- with_lecuyer_seed(2026, {
    rajiveplus:::get_random_direction_bound_robustH(
      n_obs = n,
      dims = d,
      num_samples = 30,
      num_cores = 2
    )
  })
  b <- with_lecuyer_seed(2026, {
    rajiveplus:::get_random_direction_bound_robustH(
      n_obs = n,
      dims = d,
      num_samples = 30,
      num_cores = 2
    )
  })

  expect_equal(a, b)
})

test_that("get_perm_bound_robustH is reproducible under L'Ecuyer", {
  skip_on_cran()

  set.seed(123)
  X1 <- matrix(rnorm(40 * 10), 40, 10)
  X2 <- matrix(rnorm(40 * 9), 40, 9)
  s1 <- svd(X1)
  s2 <- svd(X2)
  block_svd <- list(
    list(u = s1$u, d = s1$d, v = s1$v),
    list(u = s2$u, d = s2$d, v = s2$v)
  )
  ranks <- c(3, 3)

  a <- with_lecuyer_seed(3030, {
    rajiveplus:::get_perm_bound_robustH(
      block_svd = block_svd,
      initial_signal_ranks = ranks,
      num_samples = 25,
      num_cores = 2
    )
  })
  b <- with_lecuyer_seed(3030, {
    rajiveplus:::get_perm_bound_robustH(
      block_svd = block_svd,
      initial_signal_ranks = ranks,
      num_samples = 25,
      num_cores = 2
    )
  })

  expect_equal(a, b)
})
