context("proposal_fit = 'hybrid' with Stan fits")

# Model, data and closed-form logml: helper-hybrid-stan.R.
# The rstan version of these tests is in test-proposal-fit-hybrid-rstan.R.

test_that("CmdStanMCMC: proposal_fit = 'hybrid' uses the gradients", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  if (!file.exists(cmdstanr::cmdstan_path())) {
    testthat::skip("CmdStan is not installed in the expected path for cmdstanr.")
  }
  skip_if_rstan_loaded()

  dat <- hybrid_stan_data()
  log_true <- hybrid_stan_logml(dat)

  tf <- tempfile(fileext = ".stan")
  on.exit(unlink(tf), add = TRUE)
  writeLines(hybrid_stan_code, tf)
  mod <- cmdstanr::cmdstan_model(tf, quiet = TRUE)

  fit <- mod$sample(
    data = dat,
    seed = 909,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 1000,
    iter_sampling = 2000,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
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

  # the gradients were actually obtained and used
  testthat::expect_equal(bs_hybrid$proposal_fit_info$combiner, "geometric")
  testthat::expect_equal(
    bs_hybrid$proposal_fit_info$n_gradients_used,
    2000L
  )
  testthat::expect_equal(bs_sample$proposal_fit_info$combiner, "none")

  # and the estimate is still correct
  testthat::expect_true(is.finite(bs_hybrid$logml))
  testthat::expect_lt(abs(bs_hybrid$logml - log_true), 0.1)
  testthat::expect_lt(abs(bs_hybrid$logml - bs_sample$logml), 0.1)
})

test_that("CmdStanMCMC: the default is unchanged by the new argument", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  if (!file.exists(cmdstanr::cmdstan_path())) {
    testthat::skip("CmdStan is not installed in the expected path for cmdstanr.")
  }
  skip_if_rstan_loaded()

  dat <- hybrid_stan_data(J = 2L, N = 20L)

  tf <- tempfile(fileext = ".stan")
  on.exit(unlink(tf), add = TRUE)
  writeLines(hybrid_stan_code, tf)
  mod <- cmdstanr::cmdstan_model(tf, quiet = TRUE)
  fit <- mod$sample(
    data = dat,
    seed = 55,
    chains = 1,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )

  set.seed(1)
  a <- bridgesampling::bridge_sampler(fit, silent = TRUE, use_neff = FALSE)
  set.seed(1)
  b <- bridgesampling::bridge_sampler(
    fit,
    silent = TRUE,
    use_neff = FALSE,
    proposal_fit = "sample"
  )
  testthat::expect_equal(a$logml, b$logml)
  testthat::expect_equal(a$proposal_fit_info$combiner, "none")
})
