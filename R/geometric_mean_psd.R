## Matrix geometric mean of two positive (semi-)definite covariance
## matrices: A # B = A^{1/2} (A^{-1/2} B A^{-1/2})^{1/2} A^{1/2}.
##
## Reference: Seyboldt, Carlson, Carpenter (2026), "Preconditioning
## Hamiltonian Monte Carlo by minimizing Fisher Divergence". The
## dense Fisher-divergence-optimal preconditioner is the matrix
## geometric mean of the sample covariance of the draws and the
## inverse covariance of the scores; the helper here is the
## bridge-sampling analogue used for `proposal_fit = "hybrid"`.
##
## The geometric mean is the Riemannian midpoint of A and B on the
## positive-definite cone under the affine-invariant metric, and
## coincides with 0.5 (A + B) when A and B commute.
##
## Arguments:
##   A, B: p x p symmetric positive semi-definite matrices, already
##         nearPD-corrected by the caller.
##
## Returns: a p x p symmetric matrix. On failure (rank-deficient A,
##   non-finite output) returns the arithmetic mean 0.5 (A + B) with a
##   warning and the attribute "fallback" set to TRUE.
.geometric_mean_psd <- function(A, B) {
  if (!is.matrix(A) || !is.matrix(B))
    stop(".geometric_mean_psd(): A and B must be matrices.", call. = FALSE)
  if (!identical(dim(A), dim(B)) || nrow(A) != ncol(A))
    stop(".geometric_mean_psd(): A and B must be square and the same size.",
         call. = FALSE)

  ## Force exact symmetry up front -- both eigen() and the
  ## sandwiches below assume it.
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2

  arithmetic_fallback <- function(reason) {
    warning(sprintf(".geometric_mean_psd(): %s; falling back to ",
                    reason), "arithmetic mean.", call. = FALSE)
    out <- 0.5 * (A + B)
    out <- (out + t(out)) / 2
    attr(out, "fallback") <- TRUE
    out
  }

  eigA <- tryCatch(eigen(A, symmetric = TRUE),
                   error = function(e) NULL)
  if (is.null(eigA) || any(!is.finite(eigA$values)))
    return(arithmetic_fallback("eigendecomposition of A failed"))

  ## Clip negative eigenvalues from round-off. Small positive ones are
  ## kept: the caller pre-regularises A with nearPD(), so the geometric
  ## mean is well defined, and the non-finite check below catches any
  ## genuine numerical failure.
  lam <- eigA$values
  lam_max <- max(abs(lam))
  if (lam_max <= 0)
    return(arithmetic_fallback("A is the zero matrix"))
  lam <- pmax(lam, 0)

  U     <- eigA$vectors
  sqrtL <- sqrt(lam)
  invsL <- ifelse(sqrtL > 0, 1 / sqrtL, 0)    # pseudo-inverse for true zeros
  A_half     <- U %*% (sqrtL * t(U))
  A_neg_half <- U %*% (invsL * t(U))

  M <- A_neg_half %*% B %*% A_neg_half
  M <- (M + t(M)) / 2

  eigM <- tryCatch(eigen(M, symmetric = TRUE),
                   error = function(e) NULL)
  if (is.null(eigM) || any(!is.finite(eigM$values)))
    return(arithmetic_fallback("eigendecomposition of A^{-1/2} B A^{-1/2} failed"))

  ## Same for the inner sandwich M.
  mu <- eigM$values
  mu_max <- max(abs(mu))
  if (mu_max <= 0)
    return(arithmetic_fallback("A^{-1/2} B A^{-1/2} is the zero matrix"))
  mu <- pmax(mu, 0)

  V    <- eigM$vectors
  M_half <- V %*% (sqrt(mu) * t(V))

  G <- A_half %*% M_half %*% A_half
  G <- (G + t(G)) / 2

  if (any(!is.finite(G)))
    return(arithmetic_fallback("non-finite entries in A # B"))

  attr(G, "fallback") <- FALSE
  G
}
