
\dontrun{

################################################################################
# BAYESIAN FACTOR ANALYSIS (AS PROPOSED BY LOPES & WEST, 2004)
################################################################################

library(bridgesampling)
library(rstan)

cores <- 4
rstan_options(auto_write = TRUE)
options(mc.cores = cores)

data(ier)

#-------------------------------------------------------------------------------
# plot data
#-------------------------------------------------------------------------------

currency <- colnames(ier)
label <- c("US Dollar", "Canadian Dollar", "Yen", "Franc", "Lira", "Mark")
op <- par(mfrow = c(3, 2), mar = c(6, 6, 3, 3))

for (i in seq_along(currency)) {
  plot(ier[,currency[i]], type = "l", col = "darkblue",  axes = FALSE,
       ylim = c(-4, 4), ylab = "", xlab = "", lwd = 2)
  axis(1, at = 0:12*12, labels = 1975:1987, cex.axis = 1.7)
  axis(2, at = pretty(c(-4, 4)), las = 1, cex.axis = 1.7)
  mtext("Year", 1, cex = 1.5, line = 3.2)
  mtext("Exchange Rate Changes", 2, cex = 1.4, line = 3.2)
  mtext(label[i], 3, cex = 1.6, line = .1)
}

par(op)

#-------------------------------------------------------------------------------
# stan model
#-------------------------------------------------------------------------------

model_code <-
"data {
  int<lower=1> T; // number of observations
  int<lower=1> m; // number of variables
  int<lower=1> k; // number of factors
  matrix[T,m] Y;  // data matrix
}
transformed data {
  int<lower = 1> r;
  vector[m] zeros;
  r = m * k - k * (k - 1) / 2; // number of non-zero factor loadings
  zeros = rep_vector(0.0, m);
}
parameters {
  real beta_lower[r - k];  // lower-diagonal elements of beta
  real<lower = 0> beta_diag [k]; // diagonal elements of beta
  vector<lower = 0>[m] sigma2; // residual variances
}
transformed parameters {
  matrix[m,k] beta;
  cov_matrix[m] Omega;
  // construct lower-triangular factor loadings matrix
  {
    int index_lower = 1;
    for (j in 1:k) {
      for (i in 1:m) {
        if (i == j) {
          beta[j,j] = beta_diag[j];
        } else if (i >= j) {
          beta[i,j] = beta_lower[index_lower];
          index_lower = index_lower + 1;
        } else {
          beta[i,j] = 0.0;
        }
      }
    }
  }
  Omega = beta*beta' + diag_matrix(sigma2);
}
model {
  // priors
  target += normal_lpdf(beta_diag | 0, 1) - k * normal_lccdf(0 | 0, 1);
  target += normal_lpdf(beta_lower | 0, 1);
  target += inv_gamma_lpdf(sigma2 | 2.2 / 2.0, 0.1 / 2.0);

  // likelihood
  for(t in 1:T) {
    target += multi_normal_lpdf(Y[t] | zeros, Omega);
  }
}"

# compile model
model <- stan_model(model_code = model_code)


#-------------------------------------------------------------------------------
# fit models and compute log marginal likelihoods
#-------------------------------------------------------------------------------

# function for generating starting values
init_fun <- function(nchains, k, m) {
  r <- m * k - k * (k - 1) / 2
  out <- vector("list", nchains)
  for (i in seq_len(nchains)) {
    beta_lower <- array(runif(r - k, 0.05, 1), dim = r - k)
    beta_diag <- array(runif(k, .05, 1), dim = k)
    sigma2 <- array(runif(m, .05, 1.5), dim = m)
    out[[i]] <- list(beta_lower = beta_lower,
                     beta_diag = beta_diag,
                     sigma2 = sigma2)
  }
  return(out)
}

set.seed(1)
stanfit <- bridge <- vector("list", 3)
for (k in 1:3) {
  stanfit[[k]] <- sampling(model,
                           data = list(Y = ier, T = nrow(ier),
                                       m = ncol(ier), k = k),
                           iter = 11000, warmup = 1000, chains = 4,
                           init = init_fun(nchains = 4, k = k, m = ncol(ier)),
                           cores = 4, seed = 1)
  bridge[[k]] <- bridge_sampler(stanfit[[k]], method = "warp3",
                                repetitions = 10, cores = 4)
}

# example output
print(bridge[[2]])

#-------------------------------------------------------------------------------
# compute posterior model probabilities
#-------------------------------------------------------------------------------

post_prob(bridge[[1]], bridge[[2]], bridge[[3]],
          model_names = c("k = 1", "k = 2", "k = 3"))

}
