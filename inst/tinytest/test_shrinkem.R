# tinytest tests for shrinkem
#
# Deterministic checks (structure, reproducibility, validation) run everywhere,
# including on CRAN. Stochastic behaviour checks depend on the platform BLAS
# (different draws on Windows vs macOS vs Linux), so they are gated behind
# at_home() -- they run locally but are skipped on CRAN to avoid spurious
# cross-platform failures.

est   <- c(8, 0, -8)
Sigma <- diag(3)
it    <- 2000L
bi    <- 500L

# ===========================================================================
# Deterministic: safe on every platform (no dependence on draw values)
# ===========================================================================

set.seed(1)
fit <- shrinkem(est, Sigma, type = "horseshoe", iterations = it, burnin = bi)

expect_true(inherits(fit, "shrinkem"))
expect_equal(nrow(fit$estimates), 3L)
expect_equal(colnames(fit$estimates),
             c("input.est", "shrunk.mean", "shrunk.median", "shrunk.mode",
               "shrunk.lower", "shrunk.upper", "nonzero"))
expect_equal(dim(fit$draws$beta), c(it, 3L))            # store = 1
expect_true(all(is.finite(fit$estimates$shrunk.mean)))

# median column must not duplicate the mean (guards the old shrunk.median bug)
expect_true(any(fit$estimates$shrunk.mean != fit$estimates$shrunk.median))

# exact reproducibility on a given machine (confirms C++ uses R's RNG stream)
set.seed(42); a <- shrinkem(est, Sigma, type = "horseshoe", iterations = it, burnin = bi)
set.seed(42); b <- shrinkem(est, Sigma, type = "horseshoe", iterations = it, burnin = bi)
expect_equal(a$draws$beta, b$draws$beta)

# grouped global shrinkage: one lambda^2 column per group
set.seed(1)
fg <- shrinkem(c(8, 0, -8, 0, 8), diag(5), type = "horseshoe",
               group = c(1, 1, 2, 2, 2), iterations = it, burnin = bi)
expect_equal(ncol(fg$draws$lambda2), 2L)
expect_equal(nrow(fg$estimates), 5L)

# fixed-lambda ridge takes the i.i.d. fast path (no lambda sampling)
set.seed(1)
ff <- shrinkem(est, Sigma, type = "ridge",
               lambda2.fixed = TRUE, lambda2 = 1, iterations = it, burnin = bi)
expect_equal(dim(ff$draws$beta), c(it, 3L))
expect_true(is.null(ff$draws$lambda2))
expect_true(all(is.finite(ff$estimates$shrunk.mean)))

# all three types run and return finite estimates
set.seed(1); fl <- shrinkem(est, Sigma, type = "lasso", iterations = it, burnin = bi)
set.seed(1); fr <- shrinkem(est, Sigma, type = "ridge", iterations = it, burnin = bi)
expect_true(all(is.finite(fl$estimates$shrunk.mean)))
expect_true(all(is.finite(fr$estimates$shrunk.mean)))

# input validation
expect_error(shrinkem(est, Sigma, type = "bogus"))
expect_error(shrinkem(est, Sigma, type = "ridge", lambda2.fixed = TRUE))   # lambda2 missing
expect_error(shrinkem(est, Sigma, type = "horseshoe", group = c(1, 1)))    # wrong length

# ===========================================================================
# Stochastic behaviour: local only (skipped on CRAN via at_home())
# ===========================================================================
if (at_home()) {

  em <- fit$estimates
  # strong signals keep their sign; the null is shrunk harder than a signal
  expect_true(em$shrunk.mean[1] > 0)
  expect_true(em$shrunk.mean[3] < 0)
  expect_true(abs(em$shrunk.mean[2]) < abs(em$shrunk.mean[1]))
  expect_true(abs(em$shrunk.mode[2]) < abs(em$shrunk.mode[1]))
  # shrinkage pulls toward zero, never expands
  expect_true(all(abs(em$shrunk.mean) <= abs(em$input.est) + 0.5))
  # interval-based selection: both signals nonzero, null not
  expect_equal(em$nonzero, c(TRUE, FALSE, TRUE))

  # lasso and ridge shrink the null harder than the signal
  expect_true(abs(fl$estimates$shrunk.mean[2]) < abs(fl$estimates$shrunk.mean[1]))
  expect_true(abs(fr$estimates$shrunk.mean[2]) < abs(fr$estimates$shrunk.mean[1]))
}
