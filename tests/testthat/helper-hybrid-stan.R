# Shared by the cmdstanr and rstan files test-proposal-fit-hybrid-*.R
#
# A conjugate model whose log marginal likelihood is available in closed
# form. Each column of y is an independent group with its own mean, so
# the model has J unconstrained parameters and the marginal likelihood
# factorises over columns:
#
#   mu_j ~ N(0, tau),  y_ij | mu_j ~ N(mu_j, sigma)
#   => y_.j ~ N(0, sigma^2 I + tau^2 1 1')
#
# J > 1 matters here: with a single parameter the sample and
# score-matching covariance estimates are both scalars and the geometric
# mean has nothing to combine.
#
# The model uses target += ..._lpdf rather than ~, so that Stan keeps the
# normalizing constants and log_prob() is the log of the unnormalized
# posterior on the same scale as the closed form below.
hybrid_stan_code <- "
data {
  int<lower=1> J;
  int<lower=1> N;
  matrix[N, J] y;
  real<lower=0> sigma;
  real<lower=0> tau;
}
parameters {
  vector[J] mu;
}
model {
  target += normal_lpdf(mu | 0, tau);
  for (j in 1:J) {
    target += normal_lpdf(y[, j] | mu[j], sigma);
  }
}"

hybrid_stan_data <- function(J = 5L, N = 30L, sigma = 1, tau = 2, seed = 4321) {
  set.seed(seed)
  mu <- rnorm(J, 0, tau)
  y <- matrix(rnorm(N * J, rep(mu, each = N), sigma), nrow = N, ncol = J)
  list(J = J, N = N, y = y, sigma = sigma, tau = tau)
}

hybrid_stan_logml <- function(dat) {
  V <- dat$sigma^2 * diag(dat$N) + dat$tau^2 * matrix(1, dat$N, dat$N)
  sum(vapply(
    seq_len(dat$J),
    function(j) {
      mvtnorm::dmvnorm(dat$y[, j], mean = rep(0, dat$N), sigma = V, log = TRUE)
    },
    numeric(1)
  ))
}
