context("Iterative scheme MCSE calculations")

test_that("MCSE is finite, positive, and returned for normal method", {
  skip_if_not_installed("Brobdingnag")

  # Minimal deterministic inputs
  q11 <- c(0.10, 0.20, 0.30, 0.40)
  q12 <- c(0.02, 0.09, 0.10, 0.15)
  q21 <- c(0.05, 0.18, 0.22, 0.28)
  q22 <- c(0.01, 0.07, 0.11, 0.12)

  L <- diag(2)
  out <- bridgesampling:::.run.iterative.scheme(
    q11 = q11, q12 = q12, q21 = q21, q22 = q22,
    r0 = 1, tol = 1e-10, L = L,
    method = "normal", maxiter = 1000, silent = TRUE,
    criterion = "r", neff = length(q11), use_ess = FALSE
  )

  expect_type(out, "list")
  expect_true(is.finite(out$mcse_logml))
  expect_gt(out$mcse_logml, 0)
  # Deterministic value checks for the toy example
  expect_equal(out$logml,      0.13248979618200290, tolerance = 1e-12)
  expect_equal(out$mcse_logml, 0.02353233284705807, tolerance = 1e-12)
})

test_that("MCSE is invariant to constant shifts (warp3 vs normal)", {
  skip_if_not_installed("Brobdingnag")

  # Same base inputs as above
  q11 <- c(0.10, 0.20, 0.30, 0.40)
  q12 <- c(0.02, 0.09, 0.10, 0.15)
  q21 <- c(0.05, 0.18, 0.22, 0.28)
  q22 <- c(0.01, 0.07, 0.11, 0.12)

  L <- diag(2)

  out_normal <- bridgesampling:::.run.iterative.scheme(
    q11 = q11, q12 = q12, q21 = q21, q22 = q22,
    r0 = 1, tol = 1e-10, L = L,
    method = "normal", maxiter = 1000, silent = TRUE,
    criterion = "r", neff = length(q11), use_ess = FALSE
  )

  out_warp3 <- bridgesampling:::.run.iterative.scheme(
    q11 = q11, q12 = q12, q21 = q21, q22 = q22,
    r0 = 1, tol = 1e-10, L = L,
    method = "warp3", maxiter = 1000, silent = TRUE,
    criterion = "r", neff = length(q11), use_ess = FALSE
  )

  # warp3 adds a constant to both l1 and l2 and shifts l*, the e^(l - l*) terms are invariant
  expect_equal(out_warp3$mcse_logml, out_normal$mcse_logml, tolerance = 1e-10)
})

test_that("MCSE roughly scales like 1/sqrt(n)", {
  skip_if_not_installed("Brobdingnag")

  base_q11 <- c(0.10, 0.20, 0.30, 0.40)
  base_q12 <- c(0.02, 0.09, 0.10, 0.15)
  base_q21 <- c(0.05, 0.18, 0.22, 0.28)
  base_q22 <- c(0.01, 0.07, 0.11, 0.12)
  L <- diag(2)

  # n
  out_n <- bridgesampling:::.run.iterative.scheme(
    q11 = base_q11, q12 = base_q12, q21 = base_q21, q22 = base_q22,
    r0 = 1, tol = 1e-10, L = L,
    method = "normal", maxiter = 1000, silent = TRUE,
    criterion = "r", neff = length(base_q11), use_ess = FALSE
  )

  # 4n (replicate samples 4x)
  k <- 4
  out_4n <- bridgesampling:::.run.iterative.scheme(
    q11 = rep(base_q11, k), q12 = rep(base_q12, k),
    q21 = rep(base_q21, k), q22 = rep(base_q22, k),
    r0 = 1, tol = 1e-10, L = L,
    method = "normal", maxiter = 1000, silent = TRUE,
    criterion = "r", neff = length(base_q11) * k, use_ess = FALSE
  )

  # Expect MCSE to drop by ~ 1/sqrt(k)
  expect_lt(out_4n$mcse_logml, out_n$mcse_logml)
  expect_equal(out_4n$mcse_logml, out_n$mcse_logml / sqrt(k), tolerance = 0.05)

  # Deterministic value for the 4x replicated case
  expect_equal(out_4n$mcse_logml, 0.01052514482181162, tolerance = 1e-12)
})

test_that("Function runs with use_ess = TRUE (if posterior installed)", {
  skip_if_not_installed("Brobdingnag")
  if (!requireNamespace("posterior", quietly = TRUE)) {
    skip("posterior not installed")
  }

  q11 <- c(0.10, 0.20, 0.30, 0.40, 0.45, 0.5)
  q12 <- c(0.02, 0.09, 0.10, 0.15, 0.18, 0.2)
  q21 <- c(0.05, 0.18, 0.22, 0.28, 0.30, 0.35)
  q22 <- c(0.01, 0.07, 0.11, 0.12, 0.14, 0.15)
  L <- diag(2)

  out <- bridgesampling:::.run.iterative.scheme(
    q11 = q11, q12 = q12, q21 = q21, q22 = q22,
    r0 = 1, tol = 1e-10, L = L,
    method = "normal", maxiter = 1000, silent = TRUE,
    criterion = "r", neff = length(q11), use_ess = TRUE
  )

  expect_true(is.finite(out$mcse_logml))
  expect_gt(out$mcse_logml, 0)
})
