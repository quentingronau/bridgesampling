## ------------------------------------------------------------------------
## Example 1: Estimating the Normalizing Constant of a Two-Dimensional
##            Standard Normal Distribution
## ------------------------------------------------------------------------

library(bridgesampling)
library(mvtnorm)

samples <- rmvnorm(1e4, mean = rep(0, 2), sigma = diag(2))
colnames(samples) <- c("x1", "x2")
log_density <- function(samples.row, data) {
  -.5*t(samples.row) %*% samples.row
}

lb <- rep(-Inf, 2)
ub <- rep(Inf, 2)
names(lb) <- names(ub) <- colnames(samples)
bridge_result <- bridge_sampler(samples = samples, log_posterior = log_density,
                                data = NULL, lb = lb, ub = ub, silent = TRUE)

# compare to analytical value
analytical <- log(2*pi)
print(cbind(bridge_result$logml, analytical))

\dontrun{

## ------------------------------------------------------------------------
## Example 2: Hierarchical Normal Model
## ------------------------------------------------------------------------

# for a full description of the example, see
vignette("bridgesampling_example_jags")

library(R2jags)

### generate data ###

set.seed(12345)

mu <- 0
tau2 <- 0.5
sigma2 <- 1

n <- 20
theta <- rnorm(n, mu, sqrt(tau2))
y <- rnorm(n, theta, sqrt(sigma2))


### set prior parameters
alpha <- 1
beta <- 1
mu0 <- 0
tau20 <- 1

### functions to get posterior samples ###

### H0: mu = 0

getSamplesModelH0 <- function(data, niter = 52000, nburnin = 2000, nchains = 3) {

  model <- "
    model {
      for (i in 1:n) {
        theta[i] ~ dnorm(0, invTau2)
          y[i] ~ dnorm(theta[i], 1/sigma2)
      }
      invTau2 ~ dgamma(alpha, beta)
      tau2 <- 1/invTau2
    }"

  s <- jags(data, parameters.to.save = c("theta", "invTau2"),
            model.file = textConnection(model),
            n.chains = nchains, n.iter = niter,
            n.burnin = nburnin, n.thin = 1)

  return(s)

}

### H1: mu != 0

getSamplesModelH1 <- function(data, niter = 52000, nburnin = 2000,
                              nchains = 3) {

  model <- "
    model {
      for (i in 1:n) {
        theta[i] ~ dnorm(mu, invTau2)
        y[i] ~ dnorm(theta[i], 1/sigma2)
      }
      mu ~ dnorm(mu0, 1/tau20)
      invTau2 ~ dgamma(alpha, beta)
      tau2 <- 1/invTau2
    }"

  s <- jags(data, parameters.to.save = c("theta", "mu", "invTau2"),
            model.file = textConnection(model),
            n.chains = nchains, n.iter = niter,
            n.burnin = nburnin, n.thin = 1)

  return(s)

}

### get posterior samples ###

# create data lists for Jags
data_H0 <- list(y = y, n = length(y), alpha = alpha, beta = beta, sigma2 = sigma2)
data_H1 <- list(y = y, n = length(y), mu0 = mu0, tau20 = tau20, alpha = alpha,
                beta = beta, sigma2 = sigma2)

# fit models
samples_H0 <- getSamplesModelH0(data_H0)
samples_H1 <- getSamplesModelH1(data_H1)


### functions for evaluating the unnormalized posteriors on log scale ###
log_posterior_H0 <- function(samples.row, data) {

  mu <- 0
  invTau2 <- samples.row[[ "invTau2" ]]
  theta <- samples.row[ paste0("theta[", seq_along(data$y), "]") ]

  sum(dnorm(data$y, theta, data$sigma2, log = TRUE)) +
    sum(dnorm(theta, mu, 1/sqrt(invTau2), log = TRUE)) +
    dgamma(invTau2, data$alpha, data$beta, log = TRUE)

}

log_posterior_H1 <- function(samples.row, data) {

  mu <- samples.row[[ "mu" ]]
  invTau2 <- samples.row[[ "invTau2" ]]
  theta <- samples.row[ paste0("theta[", seq_along(data$y), "]") ]

  sum(dnorm(data$y, theta, data$sigma2, log = TRUE)) +
    sum(dnorm(theta, mu, 1/sqrt(invTau2), log = TRUE)) +
    dnorm(mu, data$mu0, sqrt(data$tau20), log = TRUE) +
    dgamma(invTau2, data$alpha, data$beta, log = TRUE)

}

# specify parameter bounds H0
cn <- colnames(samples_H0$BUGSoutput$sims.matrix)
cn <- cn[cn != "deviance"]
lb_H0 <- rep(-Inf, length(cn))
ub_H0 <- rep(Inf, length(cn))
names(lb_H0) <- names(ub_H0) <- cn
lb_H0[[ "invTau2" ]] <- 0

# specify parameter bounds H1
cn <- colnames(samples_H1$BUGSoutput$sims.matrix)
cn <- cn[cn != "deviance"]
lb_H1 <- rep(-Inf, length(cn))
ub_H1 <- rep(Inf, length(cn))
names(lb_H1) <- names(ub_H1) <- cn
lb_H1[[ "invTau2" ]] <- 0


# compute log marginal likelihood via bridge sampling for H0
H0.bridge <- bridge_sampler(samples = samples_H0, data = data_H0,
                            log_posterior = log_posterior_H0, lb = lb_H0,
                            ub = ub_H0, silent = TRUE)
print(H0.bridge)

# compute log marginal likelihood via bridge sampling for H1
H1.bridge <- bridge_sampler(samples = samples_H1, data = data_H1,
                            log_posterior = log_posterior_H1, lb = lb_H1,
                            ub = ub_H1, silent = TRUE)
print(H1.bridge)

# compute percentage error
print(error_measures(H0.bridge)$percentage)
print(error_measures(H1.bridge)$percentage)

# compute Bayes factor
BF01 <- bf(H0.bridge, H1.bridge)
print(BF01)

# compute posterior model probabilities (assuming equal prior model probabilities)
post1 <- post_prob(H0.bridge, H1.bridge)
print(post1)

# compute posterior model probabilities (using user-specified prior model probabilities)
post2 <- post_prob(H0.bridge, H1.bridge, prior_prob = c(.6, .4))
print(post2)

}

\dontrun{

## ------------------------------------------------------------------------
## Example 3: rstanarm
## ------------------------------------------------------------------------
library(rstanarm)

# N.B.: remember to specify the diagnostic_file

fit_1 <- stan_glm(mpg ~ wt + qsec + am, data = mtcars,
                  chains = 2, cores = 2, iter = 5000,
                  diagnostic_file = file.path(tempdir(), "df.csv"))
bridge_1 <- bridge_sampler(fit_1)
fit_2 <- update(fit_1, formula = . ~ . + cyl)
bridge_2 <- bridge_sampler(fit_2, method = "warp3")
bf(bridge_1, bridge_2)

}

