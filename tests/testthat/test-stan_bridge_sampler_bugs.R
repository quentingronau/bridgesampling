
context('Stan Bridge Sampler Bugs')

test_that("subscript out of bounds error", {
  ## https://github.com/quentingronau/bridgesampling/issues/26
  stan_mod = "
  data{
    int M;
    int J;
    int T;
    int E;
    int G;
    int N[G];
    int ii[M];
    int jj[M];
    int gg[M];
    int g_all[sum(N)];
    int y[M];
    matrix[J,J] obs_corr[G];
  }

  transformed data{
    int N_all = sum(N);
  }

  parameters{
    ordered[T] thresholds_raw[G,J];
    matrix<multiplier=5>[E,J] lam[G];
    matrix[N_all,E] eta;
    matrix[N_all,J] ystar_raw;
  }

  transformed parameters {
    ordered[T] thresholds[G,J];

    for(g in 1:G)
      for(j in 1:J)
        thresholds[g,j] = thresholds_raw[g,j] * 5;
  }


  model{
    matrix[N_all,J] ystar;
    int pos = 1;

    target += std_normal_lpdf(to_vector(ystar_raw));
    target += std_normal_lpdf(to_vector(eta));

    for(g in 1:G){
      int g_ids[N[g]] = segment(g_all,pos,N[g]);
      target += normal_lpdf(to_vector(eta)| 0,5);

      for(j in 1:J)
        target += std_normal_lpdf(thresholds_raw[g,j]);

      ystar[g_ids,] = eta[g_ids,] * lam[g] + ystar_raw[g_ids,];
      pos += N[g];
    }

    for(m in 1:M)
      target += ordered_logistic_lpmf(y[m] | ystar[ii[m],jj[m]],
                                             thresholds[gg[m],jj[m]]);
  }
  "
  testthat::skip_on_cran()
  testthat::skip_on_travis()
  testthat::skip_if_not_installed("rstan")
  library("rstan")
  # source("tests/testthat/test_dat.txt")
  source("test_dat.txt")

  suppressWarnings(capture.output(
    mod <- stan(model_code=stan_mod,data=test_dat, chains = 2)
  ))

  expect_warning(object = bridge_sampler(mod, silent=TRUE),
                 regexp = "Infinite value in iterative scheme, returning NA.")
})

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
    suppressWarnings(H0 <- bridge_sampler(stanfitH0, cores = 2, silent = TRUE))

    expect_s3_class(H0, "bridge")

  }
})

test_that("turtle example",{
  skip_on_cran()

  if (require(rstan)) {

    data("turtles")

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
        real alpha0 = sqrt(10.0)*alpha0_raw;
        real alpha1 = sqrt(10.0)*alpha1_raw;
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
                           data = list(y = turtles$y, x = turtles$x,
                                       nobs = nrow(turtles),
                                       m = max(turtles$clutch),
                                       clutch = turtles$clutch),
                           iter = 10500, warmup = 500, chains = 4))
  bs_m1_nc <- bridge_sampler(stanobject_m1_nc, method = "warp3",
                             repetitions = 25, silent=TRUE)

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
      real alpha0 = sqrt(10.0)*alpha0_raw;
      real alpha1 = sqrt(10.0)*alpha1_raw;
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
                           data = list(y = turtles$y, x = turtles$x,
                                       nobs = nrow(turtles),
                                       m = max(turtles$clutch),
                                       clutch = turtles$clucth),
                           iter = 10500, warmup = 500, chains = 4))

  bs_m0_nc <- bridge_sampler(stanobject_m0_nc, method = "warp3",
                             repetitions = 25, silent=TRUE)
  expect_equal(bf(bs_m0_nc, bs_m1_nc)$bf, rep(1.27, 25), tolerance = 0.02)
  }
})
