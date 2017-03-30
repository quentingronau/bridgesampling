
context('bridge_sampler.stanfit multicore works for one-parameter model.')


test_that("stan_bridge_sampler", {

  skip_on_cran()
  skip_on_travis()

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
