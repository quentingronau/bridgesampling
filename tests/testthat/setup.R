## Marginal likelihood for Stan tests
### H0: mu = 0
mH0 <- function(y, sigma2 = 1, alpha = 2, beta = 3, rel.tol = 10^(-10)) {
  n <- length(y)
  mH0integrand <- function(tau2, y, sigma2, alpha, beta) {
    (sigma2 + tau2)^(-n/2) * exp(-(n*mean(y)^2 + (n - 1)*sd(y)^2)/(2*(sigma2 + tau2))) *
      tau2^(-alpha - 1) * exp(-beta/tau2)
  }
  (2*pi)^(-n/2) * beta^alpha/gamma(alpha) * integrate(mH0integrand, 0, Inf, rel.tol = rel.tol,
                                                      y = y, sigma2 = sigma2, alpha = alpha,
                                                      beta = beta)$value
}
