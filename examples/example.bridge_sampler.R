\dontrun{

library(R2jags)
library(bridgesampling)

#########################################
###   "exact" marginal likelihoods    ###
#########################################

### H0: mu = 0

mH0 <- function(y, sigma2 = 1, alpha = 2, beta = 3, rel.tol = 10^(-40)) {

  mH0integrand <- function(tau2, y, sigma2, alpha, beta) {

    n <- length(y)

    (sigma2 + tau2)^(-n/2) * exp(-(n*mean(y)^2 + (n - 1)*sd(y)^2)/(2*(sigma2 + tau2))) *
      tau2^(-alpha - 1) * exp(-beta/tau2)

  }

  (2*pi)^(-n/2) * beta^alpha/gamma(alpha) *
    integrate(mH0integrand, 0, Inf, rel.tol = rel.tol, y = y, sigma2 = sigma2,
              alpha = alpha,beta = beta)$value

}

### H1: mu != 0

mH1 <- function(y, sigma2 = 1, mu0 = 0, tau20 = 1, alpha = 2, beta = 3,
                rel.tol = 10^(-80)) {

  mH1integrand <- function(tau2, y, sigma2, mu0, tau20, alpha, beta) {

    (sigma2 + tau2)^(-n/2) *
      exp(-1/2 * ((n*mean(y)^2 + (n - 1)*sd(y)^2)/(sigma2 + tau2) + mu0^2/tau20 -
         ((n*mean(y))/(sigma2 + tau2) +  mu0/tau20)^2 / (n/(sigma2 + tau2) + 1/tau20))) *
         (n/(sigma2 + tau2) + 1/tau20)^(-1/2) * tau2^(-alpha - 1) * exp(-beta/tau2)

  }

  (2*pi)^(-n/2) * (tau20)^(-1/2) * beta^alpha/gamma(alpha) *
    integrate(mH1integrand, 0, Inf,
              rel.tol = rel.tol, y = y,
              sigma2 = sigma2, mu0 = mu0,
              tau20 = tau20, alpha = alpha,
              beta = beta)$value

}


##########################################
### functions to get posterior samples ###
##########################################

### H0

getSamplesModelH0 <- function(y, sigma2, alpha, beta, niter = 20000,
                              nburnin = 2000, nchains = 1) {

  n <- length(y)
  data <- list("y", "n", "sigma2", "alpha", "beta")

  model <- "model {

              for (i in 1:n) {

              theta[i] ~ dnorm(0, invTau2)
              y[i] ~ dnorm(theta[i], 1/sigma2)

              }

              invTau2 ~ dgamma(alpha, beta)
              tau2 <- 1/invTau2

            }"

  s <- jags(data, parameters.to.save = c("theta", "invTau2"),
            model.file = textConnection(model),
            n.chains = nchains, n.iter = niter, n.burnin = nburnin, n.thin = 1)

  return(s)

}

### H1

getSamplesModelH1 <- function(y, sigma2, mu0, tau20, alpha, beta, niter = 20000,
                              nburnin = 2000, nchains = 1) {

  n <- length(y)
  data <- list("y", "n", "sigma2", "mu0", "tau20", "alpha", "beta")

  model <- "model {

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
            n.chains = nchains, n.iter = niter, n.burnin = nburnin, n.thin = 1)

  return(s)

}

postAsMatrix <- function(postsamples) {

  cn <- colnames(postsamples$BUGSoutput$sims.matrix)
  s <- postsamples$BUGSoutput$sims.matrix[ ,-which(cn == "deviance")] # cut off deviance
  return(s)

}


#################################################################
### functions to evaluate prior times likelihood on log scale ###
#################################################################

log.posterior.H1 <- function(samples.row, data) {

  mu <- samples.row[[ "mu" ]]
  invTau2 <- samples.row[[ "invTau2" ]]
  theta <- samples.row[ paste0("theta[", seq_along(data$y), "]") ]

  sum(dnorm(data$y, theta, data$sigma2, log = TRUE)) +
    sum(dnorm(theta, mu, 1/sqrt(invTau2), log = TRUE)) +
    dnorm(mu, data$mu0, sqrt(data$tau20), log = TRUE) +
    dgamma(invTau2, data$alpha, data$beta, log = TRUE)

}

log.posterior.H0 <- function(samples.row, data) {

  mu <- 0
  invTau2 <- samples.row[[ "invTau2" ]]
  theta <- samples.row[ paste0("theta[", seq_along(data$y), "]") ]

  sum(dnorm(data$y, theta, data$sigma2, log = TRUE)) +
    sum(dnorm(theta, mu, 1/sqrt(invTau2), log = TRUE)) +
    dgamma(invTau2, data$alpha, data$beta, log = TRUE)

}

#####################
### generate data ###
#####################

set.seed(123)

mu <- 0
tau2 <- 0.5
sigma2 <- 1

n <- 20
theta <- rnorm(n, mu, sqrt(tau2))
y <- rnorm(n, theta, sqrt(sigma2))

### prior parameters
alpha <- 2
beta <- 3
mu0 <- 0
tau20 <- 1

data <- list(y = y, mu0 = mu0, tau20 = tau20, alpha = alpha,
             beta = beta, sigma2 = sigma2)


#####################
### H0            ###
#####################

# get samples
postsamples <- getSamplesModelH0(y, sigma2 = sigma2, alpha = alpha, beta = beta,
                                 niter = 102000, nburnin = 2000, nchains = 1)
s <- postAsMatrix(postsamples)


# specify parameter bounds
lb <- rep(-Inf, ncol(s))
ub <- rep(Inf, ncol(s))
names(lb) <- names(ub) <- colnames(s)
lb[[ "invTau2" ]] <- 0

# compute log marginal likelihood via bridge sampling
H0.bridge <- bridge_sampler(samples = s, data = data,
                            log_posterior = log.posterior.H0, lb = lb,
                            ub = ub, cores = 4)

# compute "exact" marginal likelihood
mlH0.exact <- mH0(y, rel.tol = 10^(-20), sigma2 = sigma2, alpha = alpha, beta = beta)

# compare
ml.comparison.H0 <- cbind(exp(H0.bridge$logml), mlH0.exact)
colnames(ml.comparison.H0) <- c("bridge", "exact")
print(ml.comparison.H0)


#####################
### H1            ###
#####################

# get samples
postsamples <- getSamplesModelH1(y, sigma2 = sigma2, alpha = alpha, beta = beta,
                                 mu0 = mu0, tau20 = tau20, niter = 102000,
                                 nburnin = 2000, nchains = 1)
s <- postAsMatrix(postsamples)



# specify parameter bounds
lb <- rep(-Inf, ncol(s))
ub <- rep(Inf, ncol(s))
names(lb) <- names(ub) <- colnames(s)
lb[[ "invTau2" ]] <- 0

# compute log marginal likelihood via bridge sampling
H1.bridge <- bridge_sampler(samples = s, data = data,
                            log_posterior = log.posterior.H1, lb = lb,
                            ub = ub, cores = 4)

# compute "exact" marginal likelihood
mlH1.exact <- mH1(y, rel.tol = 10^(-20), sigma2 = sigma2, alpha = alpha,
                  beta = beta, mu0 = mu0, tau20 = tau20)

# compare
ml.comparison.H1 <- cbind(exp(H1.bridge$logml), mlH1.exact)
colnames(ml.comparison.H1) <- c("bridge", "exact")
print(ml.comparison.H1)


#####################
### Bayes factor  ###
#####################

# compute bridge sampling Bayes factor
BF.bridge <- exp(H0.bridge$logml - H1.bridge$logml)

# compute "exact" Bayes factor
BF.exact <- mlH0.exact/mlH1.exact

# compare
BF.comparison <- cbind(BF.bridge, BF.exact)
colnames(BF.comparison) <- c("bridge", "exact")
print(BF.comparison)

}

