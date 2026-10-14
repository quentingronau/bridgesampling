context("proposal_fit = 'hybrid' with rstan fits")

# Model, data and closed-form logml: helper-hybrid-stan.R.

test_that("stanfit: proposal_fit = 'hybrid' uses the gradients", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  # before skip_if_not_installed(): that call loads rstan
  skip_if_cmdstanr_methods_loaded()
  testthat::skip_if_not_installed("rstan")

  dat <- hybrid_stan_data()
  log_true <- hybrid_stan_logml(dat)

  sm <- suppressWarnings(
    rstan::stan_model(model_code = hybrid_stan_code, model_name = "hybrid_test")
  )
  fit <- rstan::sampling(
    sm,
    data = dat,
    seed = 909,
    chains = 2,
    iter = 3000,
    warmup = 1000,
    refresh = 0
  )

  bs_sample <- bridgesampling::bridge_sampler(
    fit,
    silent = TRUE,
    use_neff = FALSE
  )
  bs_hybrid <- bridgesampling::bridge_sampler(
    fit,
    silent = TRUE,
    use_neff = FALSE,
    proposal_fit = "hybrid"
  )

  testthat::expect_equal(bs_hybrid$proposal_fit_info$combiner, "geometric")
  testthat::expect_equal(bs_sample$proposal_fit_info$combiner, "none")
  testthat::expect_true(is.finite(bs_hybrid$logml))
  testthat::expect_lt(abs(bs_hybrid$logml - log_true), 0.1)
  testthat::expect_lt(abs(bs_hybrid$logml - bs_sample$logml), 0.1)
})
