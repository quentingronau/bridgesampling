context("matrix geometric mean of two covariance matrices")

test_that("the geometric mean of a matrix with itself is that matrix", {
  A <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  G <- .geometric_mean_psd(A, A)
  expect_equal(G, A, check.attributes = FALSE, tolerance = 1e-10)
  expect_false(attr(G, "fallback"))
})

test_that("commuting inputs give the elementwise geometric mean", {
  # Diagonal matrices commute, so A # B = diag(sqrt(a_i b_i)).
  A <- diag(c(1, 2, 3, 4))
  B <- diag(c(4, 3, 2, 1))
  expect_equal(
    .geometric_mean_psd(A, B),
    diag(sqrt(c(4, 6, 6, 4))),
    check.attributes = FALSE,
    tolerance = 1e-10
  )
})

test_that("the geometric mean is symmetric in its arguments", {
  set.seed(42)
  p <- 5
  A <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  B <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  expect_equal(
    .geometric_mean_psd(A, B),
    .geometric_mean_psd(B, A),
    check.attributes = FALSE,
    tolerance = 1e-8
  )
})

test_that("the result is symmetric and positive definite", {
  set.seed(7)
  p <- 6
  A <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  B <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  G <- .geometric_mean_psd(A, B)
  expect_equal(G, t(G), check.attributes = FALSE, tolerance = 1e-10)
  expect_gt(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("the determinant identity det(A # B) = sqrt(det A det B) holds", {
  set.seed(11)
  p <- 4
  A <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  B <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  G <- .geometric_mean_psd(A, B)
  expect_equal(
    determinant(G, logarithm = TRUE)$modulus[1],
    0.5 *
      (determinant(A, logarithm = TRUE)$modulus[1] +
        determinant(B, logarithm = TRUE)$modulus[1]),
    tolerance = 1e-8
  )
})

test_that("a degenerate input falls back to the arithmetic mean", {
  A <- matrix(0, 2, 2)
  B <- diag(2)
  expect_warning(G <- .geometric_mean_psd(A, B), "arithmetic mean")
  expect_true(attr(G, "fallback"))
  expect_equal(unclass(G), 0.5 * (A + B), check.attributes = FALSE)
})

test_that("malformed input is rejected", {
  expect_error(.geometric_mean_psd(1:4, diag(2)), "must be matrices")
  expect_error(.geometric_mean_psd(diag(2), diag(3)), "same size")
})
