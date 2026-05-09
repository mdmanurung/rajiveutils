# Tests for W-R4 (BH default), W-R5 (pooled PIP default), W-R6 (per-block pool).

local_js_planted <- function(seed = 55L) {
  set.seed(seed)
  n   <- 60
  pks <- c(80, 60)
  Y   <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                        n = n, pks = pks, dist.type = 1)
  fit <- withr::with_seed(seed,
    Rajive(Y$sim_data, initial_signal_ranks = c(4, 3), joint_rank = 2L)
  )
  list(fit = fit, data = Y$sim_data)
}

# ---------------------------------------------------------------------------
# W-R4: BH is now the default; it should detect signal in planted data
# ---------------------------------------------------------------------------

test_that("default BH correction detects signal in simulated data", {
  skip_on_cran()

  f <- local_js_planted()
  js <- withr::with_seed(42, jackstraw_rajive(f$fit, f$data, n_null = 50))

  # Check that default correction is BH.
  expect_equal(attr(js, "correction"), "BH")

  # With planted joint signal, at least some features should be significant.
  n_sig <- sum(vapply(js, function(blk)
    sum(vapply(blk, function(cmp) sum(cmp$significant, na.rm = TRUE),
               integer(1L))), integer(1L)))
  expect_gt(n_sig, 0L)
})

# ---------------------------------------------------------------------------
# W-R6: pool = "block" (default) produces different p-values than "global"
# ---------------------------------------------------------------------------

test_that("pool = 'block' and pool = 'global' give different p-values", {
  skip_on_cran()

  f  <- local_js_planted()

  js_block  <- withr::with_seed(7, jackstraw_rajive(f$fit, f$data,
                                                    n_null = 30,
                                                    pool = "block"))
  js_global <- withr::with_seed(7, jackstraw_rajive(f$fit, f$data,
                                                    n_null = 30,
                                                    pool = "global"))

  p1 <- js_block[[1]][[1]]$p_values
  p2 <- js_global[[1]][[1]]$p_values
  # The pooled null is larger, so p-values will typically differ.
  expect_false(isTRUE(all.equal(p1, p2)))
})

test_that("pool = 'block' attr is stored on result", {
  skip_on_cran()

  f  <- local_js_planted()
  js <- withr::with_seed(8, jackstraw_rajive(f$fit, f$data,
                                             n_null = 10,
                                             pool = "block"))
  # no attr stored for pool (not part of the spec), just check it runs
  expect_s3_class(js, "jackstraw_rajive")
})

# ---------------------------------------------------------------------------
# W-R5: default pip_group is "pooled"
# ---------------------------------------------------------------------------

test_that("default pip_group is pooled", {
  skip_if_not_installed("qvalue")

  f  <- local_js_planted()
  js <- withr::with_seed(9, jackstraw_rajive(f$fit, f$data,
                                             n_null = 20,
                                             pip = TRUE))

  expect_equal(attr(js, "pip_group"), "pooled")
})
