# Tests for the PIP (posterior inclusion probability) feature of
# jackstraw_rajive().
#
# Reference: Chung NC (2020). Bioinformatics, 36(10):3107-3114.
# PIP = 1 - lfdr(p-value), computed via qvalue::lfdr().

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

local_js_fixture <- function(seed = 1L, n = 60, pks = c(80, 60),
                              n_null = 20) {
  withr::with_seed(seed, {
    Y <- ajive.data.sim(K = 2, rankJ = 2, rankA = c(4, 3),
                        n = n, pks = pks, dist.type = 1)
  })
  fit <- withr::with_seed(seed,
    Rajive(Y$sim_data, initial_signal_ranks = c(4, 3))
  )
  list(fit = fit, data = Y$sim_data, n_null = n_null)
}

# ---------------------------------------------------------------------------
# Basic PIP: off by default, present when pip = TRUE
# ---------------------------------------------------------------------------

test_that("pip is absent when pip = FALSE (default)", {
  skip_if_not_installed("qvalue")

  f <- local_js_fixture()
  js <- withr::with_seed(42, jackstraw_rajive(f$fit, f$data, n_null = f$n_null))

  for (k in seq_len(attr(js, "n_blocks"))) {
    for (j in seq_len(attr(js, "joint_rank"))) {
      expect_null(js[[k]][[j]]$pip)
    }
  }
})

test_that("pip is present and numeric when pip = TRUE", {
  skip_if_not_installed("qvalue")

  f <- local_js_fixture()
  js <- withr::with_seed(42, jackstraw_rajive(f$fit, f$data,
                                              n_null = f$n_null, pip = TRUE))

  for (k in seq_len(attr(js, "n_blocks"))) {
    for (j in seq_len(attr(js, "joint_rank"))) {
      pip_vec <- js[[k]][[j]]$pip
      expect_true(is.numeric(pip_vec))
      expect_equal(length(pip_vec), length(js[[k]][[j]]$p_values))
    }
  }
})

# ---------------------------------------------------------------------------
# Range: PIPs must be in [0, 1]
# ---------------------------------------------------------------------------

test_that("PIPs are in [0, 1]", {
  skip_if_not_installed("qvalue")

  f <- local_js_fixture()
  js <- withr::with_seed(42, jackstraw_rajive(f$fit, f$data,
                                              n_null = f$n_null, pip = TRUE))

  for (k in seq_len(attr(js, "n_blocks"))) {
    for (j in seq_len(attr(js, "joint_rank"))) {
      pip_vec <- js[[k]][[j]]$pip
      expect_true(all(pip_vec >= 0 & pip_vec <= 1),
                  info = sprintf("block %d, comp %d", k, j))
    }
  }
})

# ---------------------------------------------------------------------------
# Monotonicity: lower p-value → higher PIP
# ---------------------------------------------------------------------------

test_that("PIPs are monotone decreasing in p-values (within each block/comp)", {
  skip_if_not_installed("qvalue")

  f <- local_js_fixture()
  js <- withr::with_seed(42, jackstraw_rajive(f$fit, f$data,
                                              n_null = f$n_null, pip = TRUE,
                                              pip_group = "block_component"))

  for (k in seq_len(attr(js, "n_blocks"))) {
    for (j in seq_len(attr(js, "joint_rank"))) {
      pv  <- js[[k]][[j]]$p_values
      pip <- js[[k]][[j]]$pip
      # rank correlation should be strongly negative (high p → low PIP)
      rho <- cor(pv, pip, method = "spearman")
      expect_lt(rho, 0,
                label = sprintf("Spearman rho block%d comp%d", k, j))
    }
  }
})

# ---------------------------------------------------------------------------
# pip_group: shape checks
# ---------------------------------------------------------------------------

test_that("pip_group = 'component' gives consistent shapes", {
  skip_if_not_installed("qvalue")

  f <- local_js_fixture()
  js <- withr::with_seed(42, jackstraw_rajive(f$fit, f$data,
                                              n_null = f$n_null, pip = TRUE,
                                              pip_group = "component"))

  expect_equal(attr(js, "pip_group"), "component")
  for (k in seq_len(attr(js, "n_blocks"))) {
    for (j in seq_len(attr(js, "joint_rank"))) {
      expect_equal(
        length(js[[k]][[j]]$pip),
        length(js[[k]][[j]]$p_values)
      )
    }
  }
})

test_that("pip_group = 'block_component' gives consistent shapes", {
  skip_if_not_installed("qvalue")

  f <- local_js_fixture()
  js <- withr::with_seed(42, jackstraw_rajive(f$fit, f$data,
                                              n_null = f$n_null, pip = TRUE,
                                              pip_group = "block_component"))

  expect_equal(attr(js, "pip_group"), "block_component")
  for (k in seq_len(attr(js, "n_blocks"))) {
    for (j in seq_len(attr(js, "joint_rank"))) {
      expect_equal(
        length(js[[k]][[j]]$pip),
        length(js[[k]][[j]]$p_values)
      )
    }
  }
})

test_that("pip_group = 'pooled' gives consistent shapes", {
  skip_if_not_installed("qvalue")

  f <- local_js_fixture()
  js <- withr::with_seed(42, jackstraw_rajive(f$fit, f$data,
                                              n_null = f$n_null, pip = TRUE,
                                              pip_group = "pooled"))

  expect_equal(attr(js, "pip_group"), "pooled")
  for (k in seq_len(attr(js, "n_blocks"))) {
    for (j in seq_len(attr(js, "joint_rank"))) {
      expect_equal(
        length(js[[k]][[j]]$pip),
        length(js[[k]][[j]]$p_values)
      )
    }
  }
})

# ---------------------------------------------------------------------------
# Equivalence to direct qvalue::lfdr() call (block_component grouping)
# ---------------------------------------------------------------------------

test_that("block_component PIPs match direct qvalue::lfdr() call", {
  skip_if_not_installed("qvalue")

  f <- local_js_fixture()
  js <- withr::with_seed(42, jackstraw_rajive(f$fit, f$data,
                                              n_null = f$n_null, pip = TRUE,
                                              pip_group = "block_component"))

  pv_11 <- js$block1$comp1$p_values
  expected_pip <- 1 - qvalue::lfdr(pv_11)

  expect_equal(js$block1$comp1$pip, expected_pip, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Error when qvalue missing (mocked)
# ---------------------------------------------------------------------------

test_that("pip = TRUE with missing qvalue gives a cli_abort error", {
  # Mock requireNamespace to return FALSE for qvalue
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) if (pkg == "qvalue") FALSE else TRUE,
    .package = "base"
  )
  f <- local_js_fixture()
  expect_error(
    jackstraw_rajive(f$fit, f$data, n_null = 5, pip = TRUE),
    regexp = "qvalue"
  )
})

# ---------------------------------------------------------------------------
# pip_pi0 argument is forwarded
# ---------------------------------------------------------------------------

test_that("pip_pi0 is forwarded to qvalue::lfdr()", {
  skip_if_not_installed("qvalue")

  f <- local_js_fixture()
  js_auto <- withr::with_seed(42, jackstraw_rajive(
    f$fit, f$data, n_null = f$n_null, pip = TRUE,
    pip_group = "block_component"
  ))
  js_pi0 <- withr::with_seed(42, jackstraw_rajive(
    f$fit, f$data, n_null = f$n_null, pip = TRUE,
    pip_group = "block_component", pip_pi0 = 0.8
  ))

  # With a forced pi0 != auto, results should differ
  expect_false(isTRUE(all.equal(
    js_auto$block1$comp1$pip,
    js_pi0$block1$comp1$pip
  )))
})
