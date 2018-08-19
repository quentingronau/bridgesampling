## ------------------------------------------------------------------------
library(bridgesampling)

set.seed(12345)

# Sleep data from t.test example
data(sleep)

# compute difference scores
y <- sleep$extra[sleep$group == 2] - sleep$extra[sleep$group == 1]
n <- length(y)

## ---- eval=FALSE---------------------------------------------------------
#  library(rstan)
#  
#  # models
#  stancodeH0 <- '
#  data {
#    int<lower=1> n; // number of observations
#    vector[n] y; // observations
#  }
#  parameters {
#    real<lower=0> sigma2; // variance parameter
#  }
#  model {
#    target += log(1/sigma2); // Jeffreys prior on sigma2
#    target += normal_lpdf(y | 0, sqrt(sigma2)); // likelihood
#  }
#  '
#  stancodeH1 <- '
#  data {
#    int<lower=1> n; // number of observations
#    vector[n] y; // observations
#    real<lower=0> r; // Cauchy prior scale
#  }
#  parameters {
#    real delta;
#    real<lower=0> sigma2;// variance parameter
#  }
#  model {
#    target += cauchy_lpdf(delta | 0, r); // Cauchy prior on delta
#    target += log(1/sigma2); // Jeffreys prior on sigma2
#    target += normal_lpdf(y | delta*sqrt(sigma2), sqrt(sigma2));  // likelihood
#  }
#  '
#  # compile models
#  stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
#  stanmodelH1 <- stan_model(model_code = stancodeH1, model_name="stanmodel")

## ---- eval=FALSE---------------------------------------------------------
#  # fit models
#  stanfitH0 <- sampling(stanmodelH0, data = list(y = y, n = n),
#                        iter = 20000, warmup = 1000, chains = 4, cores = 1,
#                        control = list(adapt_delta = .99))
#  stanfitH1 <- sampling(stanmodelH1, data = list(y = y, n = n, r = 1/sqrt(2)),
#                        iter = 20000, warmup = 1000, chains = 4, cores = 1,
#                        control = list(adapt_delta = .99))

## ---- echo=FALSE---------------------------------------------------------
load(system.file("extdata/", "vignette_stan_ttest.RData",
                     package = "bridgesampling"))

## ---- eval=FALSE---------------------------------------------------------
#  H0 <- bridge_sampler(stanfitH0, silent = TRUE)
#  H1 <- bridge_sampler(stanfitH1, silent = TRUE)

## ------------------------------------------------------------------------
print(H0)
print(H1)

## ----eval=FALSE----------------------------------------------------------
#  # compute percentage errors
#  H0.error <- error_measures(H0)$percentage
#  H1.error <- error_measures(H1)$percentage

## ------------------------------------------------------------------------
print(H0.error)
print(H1.error)

## ------------------------------------------------------------------------
# compute Bayes factor
BF10 <- bf(H1, H0)
print(BF10)

## ---- message=FALSE------------------------------------------------------
library(BayesFactor)
print(ttestBF(y))

## ---- eval=FALSE---------------------------------------------------------
#  stancodeHplus <- '
#  data {
#    int<lower=1> n; // number of observations
#    vector[n] y; // observations
#    real<lower=0> r; // Cauchy prior scale
#  }
#  parameters {
#    real<lower=0> delta; // constrained to be positive
#    real<lower=0> sigma2;// variance parameter
#  }
#  model {
#    target += cauchy_lpdf(delta | 0, r) - cauchy_lccdf(0 | 0, r); // Cauchy prior on delta
#    target += log(1/sigma2); // Jeffreys prior on sigma2
#    target += normal_lpdf(y | delta*sqrt(sigma2), sqrt(sigma2));  // likelihood
#  }
#  '
#  # compile and fit model
#  stanmodelHplus <- stan_model(model_code = stancodeHplus, model_name="stanmodel")
#  stanfitHplus <- sampling(stanmodelHplus, data = list(y = y, n = n, r = 1/sqrt(2)),
#                           iter = 30000, warmup = 1000, chains = 4,
#                           control = list(adapt_delta = .99))

## ----eval=FALSE----------------------------------------------------------
#  Hplus <- bridge_sampler(stanfitHplus, silent = TRUE)

## ------------------------------------------------------------------------
print(Hplus)

## ----eval=FALSE----------------------------------------------------------
#  Hplus.error <- error_measures(Hplus)$percentage

## ------------------------------------------------------------------------
print(Hplus.error)

## ------------------------------------------------------------------------
# compute Bayes factor
BFplus0 <- bf(Hplus, H0)
print(BFplus0)

## ------------------------------------------------------------------------
print(ttestBF(y, nullInterval = c(0, Inf)))

