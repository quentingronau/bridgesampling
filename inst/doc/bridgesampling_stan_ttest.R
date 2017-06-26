## ------------------------------------------------------------------------
library(bridgesampling)

# Sleep data from t.test example
data(sleep)

# compute difference scores
y <- sleep$extra[sleep$group == 2] - sleep$extra[sleep$group == 1]
n <- length(y)

## ---- message = FALSE, results='hide'------------------------------------
library(rstan)

# models
stancodeH0 <- '
data {
  int<lower=1> n; // number of observations
  vector[n] y; // observations
}
parameters {
  real<lower=0> sigma2; // variance parameter
}
model {
  target += log(1/sigma2); // Jeffreys prior on sigma2
  target += normal_lpdf(y | 0, sqrt(sigma2)); // likelihood
}
'
stancodeH1 <- '
data {
  int<lower=1> n; // number of observations
  vector[n] y; // observations
  real<lower=0> r; // Cauchy prior scale
}
parameters {
  real delta;
  real<lower=0> sigma2;// variance parameter
}
model {
  target += cauchy_lpdf(delta | 0, r); // Cauchy prior on delta
  target += log(1/sigma2); // Jeffreys prior on sigma2
  target += normal_lpdf(y | delta*sqrt(sigma2), sqrt(sigma2));  // likelihood
}
'
# compile models
stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
stanmodelH1 <- stan_model(model_code = stancodeH1, model_name="stanmodel")

## ---- message = FALSE, results='hide', warning=FALSE---------------------
# fit models
stanfitH0 <- sampling(stanmodelH0, data = list(y = y, n = n),
                      iter = 20000, warmup = 1000, chains = 4, cores = 1,
                         control = list(adapt_delta = .99))
stanfitH1 <- sampling(stanmodelH1, data = list(y = y, n = n, r = 1/sqrt(2)),
                      iter = 20000, warmup = 1000, chains = 4, cores = 1,
                         control = list(adapt_delta = .99))

## ------------------------------------------------------------------------
set.seed(12345)
H0 <- bridge_sampler(stanfitH0, silent = TRUE)
print(H0)
H1 <- bridge_sampler(stanfitH1, silent = TRUE)
print(H1)

## ------------------------------------------------------------------------
# compute percentage errors
print(error_measures(H0)$percentage)
print(error_measures(H1)$percentage)

## ------------------------------------------------------------------------
# compute Bayes factor
BF10 <- bf(H1, H0)
print(BF10)

## ---- message=FALSE------------------------------------------------------
library(BayesFactor)
print(ttestBF(y))

## ---- message = FALSE, results='hide', warnings = FALSE------------------
stancodeHplus <- '
data {
  int<lower=1> n; // number of observations
  vector[n] y; // observations
  real<lower=0> r; // Cauchy prior scale
}
parameters {
  real<lower=0> delta; // constrained to be positive
  real<lower=0> sigma2;// variance parameter
}
model {
  target += cauchy_lpdf(delta | 0, r) - cauchy_lccdf(0 | 0, r); // Cauchy prior on delta
  target += log(1/sigma2); // Jeffreys prior on sigma2
  target += normal_lpdf(y | delta*sqrt(sigma2), sqrt(sigma2));  // likelihood
}
'
# compile models
stanmodelHplus <- stan_model(model_code = stancodeHplus, model_name="stanmodel")
stanfitHplus <- sampling(stanmodelHplus, data = list(y = y, n = n, r = 1/sqrt(2)),
                         iter = 30000, warmup = 1000, chains = 4,
                         control = list(adapt_delta = .99))

## ------------------------------------------------------------------------
Hplus <- bridge_sampler(stanfitHplus, silent = TRUE)
print(Hplus)
print(error_measures(Hplus)$percentage)

## ------------------------------------------------------------------------
# compute Bayes factor
BFplus0 <- bf(Hplus, H0)
print(BFplus0)

## ------------------------------------------------------------------------
print(ttestBF(y, nullInterval = c(0, Inf)))

