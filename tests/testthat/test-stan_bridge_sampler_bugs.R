
context('Stan Bridge Sampler Bugs')


test_that("bridge_sampler.stanfit multicore works for one-parameter model.", {

  skip_on_cran()
  skip_on_travis()
  skip_on_os("windows")

  if (require(rstan)) {
    set.seed(12345)

    # compute difference scores
    n <- 10
    y <- rnorm(n)

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
    # compile models
    tmp <- capture.output(
    stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
    )
    # fit models
    tmp <- capture.output(
    stanfitH0 <- sampling(stanmodelH0, data = list(y = y, n = n),
                          iter = 10000, warmup = 1000, chains = 4,
                          control = list(adapt_delta = 0.95))
    )
    ######### bridge sampling ###########
    H0 <- bridge_sampler(stanfitH0, cores = 2, silent = TRUE)

  }
})

test_that("turtle example",{
  skip_on_cran()

  if (require(rstan)) {
    #-------------------------------------------------------------------------------
    # turtle data (obtained from Overstall & Forster, 2010)
    #-------------------------------------------------------------------------------

    clutch <- c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,
                5,5,5,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,
                9,9,9,9,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,
                12,12,13,13,13,14,15,15,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,
                18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,19,20,
                20,20,20,20,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,23,
                23,23,23,23,23,23,23,23,23,23,23,24,25,25,25,25,25,25,25,25,25,25,26,
                26,26,26,26,26,26,26,27,27,27,27,27,27,27,28,28,28,28,28,28,28,28,29,
                29,29,29,29,29,30,30,30,30,30,30,30,30,30,30,30,30,30,30,31,31,31,31,
                31,31,31,31,31,31,31,31,31,31,31,31,31,31)

    bwt <- c(5.66, 5.66, 5.45, 5.49, 5.70, 6.06, 6.23, 5.64, 5.28, 5.59, 6.15,
             6.23, 5.96, 6.51, 5.73, 5.20, 6.71, 5.55, 5.55, 6.39, 6.46, 6.41,
             6.81, 6.84, 6.36, 6.42, 5.98, 6.52, 6.78, 4.78, 4.64, 4.95, 5.69,
             4.87, 4.13, 5.04, 4.81, 8.00, 6.36, 7.95, 7.76, 7.64, 7.72, 6.77,
             7.90, 6.84, 7.30, 7.85, 7.47, 7.00, 6.84, 7.61, 7.29, 7.32, 8.66,
             7.96, 7.63, 7.01, 7.77, 8.49, 7.11, 7.76, 8.79, 8.14, 7.63, 5.58,
             7.21, 7.44, 6.98, 6.05, 8.36, 8.53, 6.32, 6.02, 6.50, 6.68, 6.47,
             5.55, 6.74, 6.31, 6.26, 4.90, 5.46, 5.21, 6.70, 6.36, 5.99, 4.98,
             6.64, 5.90, 5.61, 6.30, 7.10, 6.68, 8.10, 5.41, 5.11, 7.46, 5.83,
             5.57, 3.83, 4.15, 5.97, 7.36, 7.66, 4.83, 8.19, 7.60, 7.55, 8.37,
             5.86, 7.69, 6.33, 6.67, 7.45, 6.44, 4.33, 7.29, 6.22, 6.05, 4.22,
             4.91, 8.00, 5.96, 7.55, 6.31, 6.44, 6.49, 5.54, 5.93, 7.50, 6.10,
             6.22, 6.81, 6.48, 6.00, 5.92, 4.59, 4.85, 5.47, 4.46, 6.21, 7.20,
             8.08, 7.81, 8.11, 7.73, 7.92, 5.55, 8.07, 8.02, 4.79, 4.72, 6.00,
             7.36, 5.61, 5.70, 5.89, 5.75, 4.09, 3.18, 6.14, 6.39, 4.14, 5.08,
             6.04, 5.93, 5.53, 4.89, 4.32, 6.75, 6.49, 7.67, 5.25, 5.41, 5.15,
             5.28, 5.85, 6.79, 6.17, 5.25, 5.06, 5.62, 4.92, 6.29, 6.33, 6.13,
             6.11, 5.88, 5.28, 5.66, 7.45, 6.93, 6.99, 7.56, 7.70, 8.48, 7.87,
             4.97, 5.83, 7.07, 6.63, 6.43, 6.30, 6.83, 5.70, 5.99, 6.34, 6.01,
             6.47, 7.29, 6.02, 6.10, 5.09, 6.38, 5.15, 6.00, 5.96, 5.66, 6.05,
             5.71, 5.18, 6.65, 6.28, 6.56, 5.11, 5.99, 6.30, 6.00, 6.35, 5.40,
             5.39, 5.50, 6.52, 7.88, 7.15, 6.23, 6.71, 5.87, 6.51, 6.97, 6.10,
             5.58, 5.20)

    surv <- c(1,1,0,0,0,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
              0,1,1,0,1,1,0,1,0,1,1,0,1,0,1,0,0,0,0,1,1,1,0,0,1,0,0,1,1,1,1,1,0,1,1,
              1,1,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,1,0,1,1,0,1,1,1,0,0,0,1,1,0,0,0,0,1,
              1,0,0,0,1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,
              0,1,1,1,1,1,0,1,0,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,
              0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,0,1,1,0,1,0,0,0,0,0,1,0,1,0,1,
              1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0)

    nobs <- length(surv)
    m <- max(clutch)
    # sort clutches by mean birthweight
    means <- tapply(bwt, clutch, mean)
    index <- sort(means, index.return = TRUE)$ix
    index_2 <- sapply(unique(clutch), function(x) which(index == x))
    clutch_o <- index_2[clutch]
    ### m1 (model with random intercepts) ###
    m1_code_nc <-
      "data {
        int<lower = 1> nobs;
        int<lower = 0, upper = 1> y[nobs];
        real<lower = 0> x[nobs];
        int<lower = 1> m;
        int<lower = 1> clutch[nobs];
      }
      parameters {
        real alpha0_raw;
        real alpha1_raw;
        vector[m] b_raw;
        real<lower = 0> sigma2;
      }
      transformed parameters {
        vector[m] b;
        real<lower = 0> sigma = sqrt(sigma2);
        real alpha0 = sqrt(10)*alpha0_raw;
        real alpha1 = sqrt(10)*alpha1_raw;
        b = b_raw*sigma;
      }
      model {
        // priors
        target += -2*log(1 + sigma2); // p(sigma2) = 1/(1 + sigma2)^2
        target += normal_lpdf(alpha0_raw | 0, 1);
        target += normal_lpdf(alpha1_raw | 0, 1);

        // random effects
        target += normal_lpdf(b_raw | 0, 1);

        // likelihood
        for (i in 1:nobs)
        target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1*x[i] + b[clutch[i]]));
  }"

  tmp <- capture.output(stanobject_m1_nc <- stan(model_code = m1_code_nc,
                           data = list(y = surv, x = bwt, nobs = nobs,
                                       m = m, clutch = clutch_o),
                           iter = 10500, warmup = 500, chains = 4))
  bs_m1_nc <- bridge_sampler(stanobject_m1_nc, method = "warp3", repetitions = 25, silent=TRUE)

  m0_code_nc <-
    "data {
      int<lower = 1> nobs;
      int<lower = 0, upper = 1> y[nobs];
      real<lower = 0> x[nobs];
    }
    parameters {
      real alpha0_raw;
      real alpha1_raw;
    }
    transformed parameters {
      real alpha0 = sqrt(10)*alpha0_raw;
      real alpha1 = sqrt(10)*alpha1_raw;
    }
    model {
      // priors
      target += normal_lpdf(alpha0_raw | 0, 1);
      target += normal_lpdf(alpha1_raw | 0, 1);

      // likelihood
      for (i in 1:nobs)
        target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1*x[i]));
    }"

  tmp <- capture.output(stanobject_m0_nc <- stan(model_code = m0_code_nc,
                           data = list(y = surv, x = bwt, nobs = nobs,
                                       m = m, clutch = clutch_o),
                           iter = 10500, warmup = 500, chains = 4))

  bs_m0_nc <- bridge_sampler(stanobject_m0_nc, method = "warp3", repetitions = 25, silent=TRUE)
  expect_equal(bf(bs_m0_nc, bs_m1_nc)$bf, rep(1.27, 25), tolerance = 0.02)
  }
})
