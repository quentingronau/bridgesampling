## ------------------------------------------------------------------------
library(bridgesampling)

### generate data ###
set.seed(12345)

mu <- 0
tau2 <- 0.5
sigma2 <- 1

n <- 20
theta <- rnorm(n, mu, sqrt(tau2))
y <- rnorm(n, theta, sqrt(sigma2))
  

## ----eval=FALSE----------------------------------------------------------
#  ### set prior parameters ###
#  mu0 <- 0
#  tau20 <- 1
#  alpha <- 1
#  beta <- 1

## ---- eval=FALSE---------------------------------------------------------
#  library(rstan)
#  
#  # models
#  stancodeH0 <- 'data {
#    int<lower=1> n; // number of observations
#    vector[n] y; // observations
#    real<lower=0> alpha;
#    real<lower=0> beta;
#    real<lower=0> sigma2;
#  }
#  parameters {
#    real<lower=0> tau2; // group-level variance
#    vector[n] theta; // participant effects
#  }
#  model {
#    target += inv_gamma_lpdf(tau2 | alpha, beta);
#    target += normal_lpdf(theta | 0, sqrt(tau2));
#    target += normal_lpdf(y | theta, sqrt(sigma2));
#  }
#  '
#  stancodeH1 <- 'data {
#    int<lower=1> n; // number of observations
#    vector[n] y; // observations
#    real mu0;
#    real<lower=0> tau20;
#    real<lower=0> alpha;
#    real<lower=0> beta;
#    real<lower=0> sigma2;
#  }
#  parameters {
#    real mu;
#    real<lower=0> tau2; // group-level variance
#    vector[n] theta; // participant effects
#  }
#  model {
#    target += normal_lpdf(mu | mu0, sqrt(tau20));
#    target += inv_gamma_lpdf(tau2 | alpha, beta);
#    target += normal_lpdf(theta | mu, sqrt(tau2));
#    target += normal_lpdf(y | theta, sqrt(sigma2));
#  }
#  '
#  # compile models
#  stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
#  stanmodelH1 <- stan_model(model_code = stancodeH1, model_name="stanmodel")

## ---- eval=FALSE---------------------------------------------------------
#  # fit models
#  stanfitH0 <- sampling(stanmodelH0, data = list(y = y, n = n,
#                                                 alpha = alpha,
#                                                 beta = beta,
#                                                 sigma2 = sigma2),
#                        iter = 50000, warmup = 1000, chains = 3, cores = 1)
#  stanfitH1 <- sampling(stanmodelH1, data = list(y = y, n = n,
#                                                 mu0 = mu0,
#                                                 tau20 = tau20,
#                                                 alpha = alpha,
#                                                 beta = beta,
#                                                 sigma2 = sigma2),
#                        iter = 50000, warmup = 1000, chains = 3, cores = 1)

## ---- echo=FALSE---------------------------------------------------------
load(system.file("extdata/", "vignette_example_stan.RData",
                     package = "bridgesampling"))

## ----eval=FALSE----------------------------------------------------------
#  # compute log marginal likelihood via bridge sampling for H0
#  H0.bridge <- bridge_sampler(stanfitH0, silent = TRUE)
#  
#  # compute log marginal likelihood via bridge sampling for H1
#  H1.bridge <- bridge_sampler(stanfitH1, silent = TRUE)

## ------------------------------------------------------------------------
print(H0.bridge)
print(H1.bridge)

## ----eval=FALSE----------------------------------------------------------
#  # compute percentage errors
#  H0.error <- error_measures(H0.bridge)$percentage
#  H1.error <- error_measures(H1.bridge)$percentage

## ------------------------------------------------------------------------
print(H0.error)
print(H1.error)

## ------------------------------------------------------------------------
# compute Bayes factor
BF01 <- bf(H0.bridge, H1.bridge)
print(BF01)

## ------------------------------------------------------------------------
# compute posterior model probabilities (assuming equal prior model probabilities)
post1 <- post_prob(H0.bridge, H1.bridge)
print(post1)

## ------------------------------------------------------------------------
# compute posterior model probabilities (using user-specified prior model probabilities)
post2 <- post_prob(H0.bridge, H1.bridge, prior_prob = c(.6, .4))
print(post2)

