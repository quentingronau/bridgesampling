
\dontrun{

library(bridgesampling)
library(rstan)

data(turtles)

#-------------------------------------------------------------------------------
# plot data
#-------------------------------------------------------------------------------

# reproduce Figure 1 from Sinharay & Stern (2005)
xticks <- pretty(turtles$clutch)
yticks <- pretty(turtles$x)

plot(1, type = "n", axes = FALSE, ylab = "", xlab = "", xlim = range(xticks),
     ylim =  range(yticks))
points(turtles$clutch, turtles$x, pch = ifelse(turtles$y, 21, 4), cex = 1.3,
       col = ifelse(turtles$y, "black", "darkred"), bg = "grey", lwd = 1.3)
axis(1, cex.axis = 1.4)
mtext("Clutch Identifier", side = 1, line = 2.9, cex = 1.8)
axis(2, las = 1, cex.axis = 1.4)
mtext("Birth Weight (Grams)", side = 2, line = 2.6, cex = 1.8)

#-------------------------------------------------------------------------------
# Analysis: Natural Selection Study (compute same BF as Sinharay & Stern, 2005)
#-------------------------------------------------------------------------------

### H0 (model without random intercepts) ###
H0_code <-
"data {
  int<lower = 1> N;
  int<lower = 0, upper = 1> y[N];
  real<lower = 0> x[N];
}
parameters {
  real alpha0_raw;
  real alpha1_raw;
}
transformed parameters {
  real alpha0 = sqrt(10.0) * alpha0_raw;
  real alpha1 = sqrt(10.0) * alpha1_raw;
}
model {
  // priors
  target += normal_lpdf(alpha0_raw | 0, 1);
  target += normal_lpdf(alpha1_raw | 0, 1);

  // likelihood
  for (i in 1:N) {
    target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1 * x[i]));
  }
}"

### H1 (model with random intercepts) ###
H1_code <-
"data {
  int<lower = 1> N;
  int<lower = 0, upper = 1> y[N];
  real<lower = 0> x[N];
  int<lower = 1> C;
  int<lower = 1, upper = C> clutch[N];
}
parameters {
  real alpha0_raw;
  real alpha1_raw;
  vector[C] b_raw;
  real<lower = 0> sigma2;
}
transformed parameters {
  vector[C] b;
  real<lower = 0> sigma = sqrt(sigma2);
  real alpha0 = sqrt(10.0) * alpha0_raw;
  real alpha1 = sqrt(10.0) * alpha1_raw;
  b = sigma * b_raw;
}
model {
  // priors
  target += - 2 * log(1 + sigma2); // p(sigma2) = 1 / (1 + sigma2) ^ 2
  target += normal_lpdf(alpha0_raw | 0, 1);
  target += normal_lpdf(alpha1_raw | 0, 1);

  // random effects
  target += normal_lpdf(b_raw | 0, 1);

  // likelihood
  for (i in 1:N) {
    target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1 * x[i] + b[clutch[i]]));
  }
}"

set.seed(1)
### fit models ###
stanfit_H0 <- stan(model_code = H0_code,
                   data = list(y = turtles$y, x = turtles$x, N = nrow(turtles)),
                   iter = 15500, warmup = 500, chains = 4, seed = 1)
stanfit_H1 <- stan(model_code = H1_code,
                   data = list(y = turtles$y, x = turtles$x, N = nrow(turtles),
                               C = max(turtles$clutch), clutch = turtles$clutch),
                   iter = 15500, warmup = 500, chains = 4, seed = 1)

set.seed(1)
### compute (log) marginal likelihoods ###
bridge_H0 <- bridge_sampler(stanfit_H0)
bridge_H1 <- bridge_sampler(stanfit_H1)

### compute approximate percentage errors ###
error_measures(bridge_H0)$percentage
error_measures(bridge_H1)$percentage

### compute Bayes factor ("true" value: BF01 = 1.273) ###
bf(bridge_H0, bridge_H1)

}
