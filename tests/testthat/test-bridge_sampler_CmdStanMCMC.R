context("test bridge_sampler cmdstanmcmc method")

testthat::test_that("bridge_sampler() works for CmdStanMCMC and basic sanity checks", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("bridgesampling")
  testthat::skip_if_not_installed("posterior")

  # Require a working CmdStan toolchain
  if (!file.exists(cmdstanr::cmdstan_path())) {
    testthat::skip("CmdStan is not installed in the expected path for cmdstanr.")
  }

  set.seed(123)

  N     <- 60L
  sigma <- 1
  y     <- rnorm(N, mean = 0.5, sd = sigma)

  data_list <- list(N = N, y = y, sigma = sigma)

  stan_code <- "
  data {
    int<lower=1> N;
    vector[N] y;
    real<lower=0> sigma;
  }
  parameters {
    real mu;
  }
  model {
    mu ~ normal(0, 1);
    y ~ normal(mu, sigma);
  }"

  tf <- tempfile(fileext = ".stan")
  on.exit(unlink(tf), add = TRUE)
  writeLines(stan_code, tf)
  mod <- cmdstanr::cmdstan_model(tf, quiet = TRUE, force_recompile = TRUE)

  fit <- mod$sample(
    data = data_list,
    seed = 202,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 2000,
    iter_sampling = 4000,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )

  bs <- bridgesampling::bridge_sampler(fit, silent = TRUE, use_neff = FALSE)

  testthat::expect_s3_class(fit, "CmdStanMCMC")
  testthat::expect_true(is.list(bs))
  testthat::expect_true(is.finite(bs$logml))

  fit2 <- mod$sample(
    data = data_list,
    seed = 203,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 750,
    iter_sampling = 2000,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
  bs2 <- bridgesampling::bridge_sampler(fit2, silent = TRUE, use_neff = FALSE)
  testthat::expect_true(is.finite(bs2$logml))
})

testthat::test_that("CmdStanMCMC bridge estimate roughly agrees with rstan", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("bridgesampling")
  if (!requireNamespace("rstan", quietly = TRUE)) testthat::skip("rstan not installed")

  if (!file.exists(cmdstanr::cmdstan_path())) {
    testthat::skip("CmdStan is not installed in the expected path for cmdstanr.")
  }

  set.seed(456)

  # Same model/data as above
  N     <- 60L
  sigma <- 1
  y     <- rnorm(N, mean = 0.25, sd = sigma)
  data_list <- list(N = N, y = y, sigma = sigma)

  stan_code <- "
  data {
    int<lower=1> N;
    vector[N] y;
    real<lower=0> sigma;
  }
  parameters {
    real mu;
  }
  model {
    mu ~ normal(0, 1);
    y ~ normal(mu, sigma);
  }"

  # --- CmdStanR/RStan fit ---
  tf <- tempfile(fileext = ".stan")
  on.exit(unlink(tf), add = TRUE)
  writeLines(stan_code, tf)
  mod_cs <- cmdstanr::cmdstan_model(tf, quiet = TRUE)
  fit_cs <- mod_cs$sample(
    data = data_list,
    seed = 777,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 10000,
    iter_sampling = 3000,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
  bs_cmd <- bridgesampling::bridge_sampler(fit_cs, silent = TRUE, use_neff = FALSE)
  testthat::expect_true(is.finite(bs_cmd$logml))

  sm <- rstan::stan_model(model_code = stan_code)
  fit_rs <- rstan::sampling(
    sm, data = data_list, seed = 777,
    chains = 4, iter = 10000, warmup = 3000, refresh = 0
  )
  bs_rstan <- bridgesampling::bridge_sampler(fit_rs, silent = TRUE, use_neff = FALSE)
  testthat::expect_true(is.finite(bs_rstan$logml))

  # Compare the two bridge estimates. Tolerance accounts for MC/bridge variance.
  testthat::expect_lt(abs(bs_cmd$logml - bs_rstan$logml), 0.5)
})
