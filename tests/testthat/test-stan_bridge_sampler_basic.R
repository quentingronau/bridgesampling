
context('stan_bridge_sampler works.')

test_that("stan_bridge_sampler", {
  if (require(rstan)) {
    set.seed(12345)

    mu <- 0
    tau2 <- 0.5
    sigma2 <- 1

    n <- 20
    theta <- rnorm(n, mu, sqrt(tau2))
    y <- rnorm(n, theta, sqrt(sigma2))

    ### set prior parameters ###
    mu0 <- 0
    tau20 <- 1
    alpha <- 1
    beta <- 1

    # models
    stancodeH0 <- 'data {
    int<lower=1> n; // number of observations
    vector[n] y; // observations
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0> sigma2;
    }
    parameters {
    real<lower=0> tau2; // group-level variance
    vector[n] theta; // participant effects
    }
    model {
    target += inv_gamma_lpdf(tau2 | alpha, beta);
    target += normal_lpdf(theta | 0, sqrt(tau2));
    target += normal_lpdf(y | theta, sqrt(sigma2));
    }
    '

    # compile models
    tmp <- capture.output(suppressMessages(
    stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
    ))

    # fit models
    tmp <- capture.output(
    stanobjectH0 <- sampling(stanmodelH0, data = list(y = y, n = n,
                                                      alpha = alpha,
                                                      beta = beta),
                             iter = 2500, warmup = 500, chains = 4, show_messages = FALSE))
    expect_is(
    H0_bridge_norm <- stan_bridge_sampler(stanobjectH0, method = "normal", silent = TRUE)
    , "bridge")

    expect_is(
    H0_bridge_warp3 <- stan_bridge_sampler(stanobjectH0, method = "warp3", silent = TRUE)
    , "bridge")

  }
})


test_that("stan_bridge_sampler in multicore", {
  testthat::skip_on_cran()
  testthat::skip_on_travis()
  if (require(rstan)) {
    set.seed(12345)

    mu <- 0
    tau2 <- 0.5
    sigma2 <- 1

    n <- 20
    theta <- rnorm(n, mu, sqrt(tau2))
    y <- rnorm(n, theta, sqrt(sigma2))

    ### set prior parameters ###
    mu0 <- 0
    tau20 <- 1
    alpha <- 1
    beta <- 1

    # models
    stancodeH0 <- 'data {
    int<lower=1> n; // number of observations
    vector[n] y; // observations
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0> sigma2;
    }
    parameters {
    real<lower=0> tau2; // group-level variance
    vector[n] theta; // participant effects
    }
    model {
    target += inv_gamma_lpdf(tau2 | alpha, beta);
    target += normal_lpdf(theta | 0, sqrt(tau2));
    target += normal_lpdf(y | theta, sqrt(sigma2));
    }
    '

    # compile models
    tmp <- capture.output(suppressMessages(
    stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
    ))

    # fit models
    tmp <- capture.output(
    stanobjectH0 <- sampling(stanmodelH0, data = list(y = y, n = n,
                                                      alpha = alpha,
                                                      beta = beta),
                             iter = 2500, warmup = 500, chains = 4, show_messages = FALSE))
    expect_is(
    H0_bridge_norm <- stan_bridge_sampler(stanobjectH0, method = "normal", silent = TRUE, cores = 2)
    , "bridge")

    expect_is(
    H0_bridge_warp3 <- stan_bridge_sampler(stanobjectH0, method = "warp3", silent = TRUE, cores = 2)
    , "bridge")

  }
})
