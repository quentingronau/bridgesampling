context(".cmdstan_log_posterior helper function for CmdStanMCMC objects")

testthat::test_that(".cmdstan_log_posterior matches lp__ and is wired the same way bridge_sampler uses it", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("posterior")
  testthat::skip_if_not_installed("bridgesampling")

  if (!file.exists(cmdstanr::cmdstan_path())) {
    testthat::skip("CmdStan is not installed in the expected path for cmdstanr.")
  }

  pkg <- "bridgesampling"
  testthat::expect_true(
    exists(".cmdstan_log_posterior", envir = asNamespace(pkg), inherits = FALSE),
    info = "Internal function .cmdstan_log_posterior not found in 'bridgesampling'."
  )
  .cmdstan_log_posterior <- get(".cmdstan_log_posterior", envir = asNamespace(pkg))

  set.seed(321)
  N     <- 40L
  sigma <- 1
  y     <- rnorm(N, 0.25, sigma)
  data_list <- list(N = N, y = y, sigma = sigma)

  stan_code <- "
  data { int<lower=1> N; vector[N] y; real<lower=0> sigma; }
  parameters { real mu; }
  model { mu ~ normal(0, 1); y ~ normal(mu, sigma); }"

  tf <- withr::local_tempfile(fileext = ".stan")
  writeLines(stan_code, tf)

  mod <- cmdstanr::cmdstan_model(tf, quiet = TRUE, force_recompile = TRUE)

  fit <- mod$sample(
    data = data_list,
    seed = 404,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 300,
    iter_sampling = 800,
    refresh = 0
  )

  # CmdStan's lp__
  draws_df <- fit$draws(variables = "lp__", format = "df")
  testthat::expect_true(nrow(draws_df) > 0)
  lp_vec_all <- draws_df$lp__
  n_all <- length(lp_vec_all)
  testthat::expect_true(n_all > 10)

  # How bridge_sampler splits internally 
  n_fit <- round(n_all / 2)
  idx_iter <- seq.int(n_fit + 1L, n_all)  # indices corresponding to q11

  # Run the actual method that wires the helper into the pipeline
  bs_out <- bridgesampling::bridge_sampler(
    samples     = fit,
    repetitions = 1,
    method      = "normal",
    cores       = 1L,
    use_neff    = FALSE,  
    silent      = TRUE,
    verbose     = FALSE
  )

  # q11 are the posterior log-densities for the iterative half, should match lp__
  testthat::expect_equal(length(bs_out$q11), length(idx_iter))
  testthat::expect_equal(
    unname(bs_out$q11),
    unname(lp_vec_all[idx_iter]),
    tolerance = 1e-4
  )

  # The helper is applied row-wise to unconstrained parameters, with `data = fit`.
  upars <- fit$unconstrain_draws(format = "matrix")
  testthat::expect_true(nrow(upars) == n_all)

  # pick a few rows from the iterative half to validate direct calls
  take <- idx_iter[seq_len(min(5L, length(idx_iter)))]
  direct_vals <- apply(upars[take, , drop = FALSE], 1, .cmdstan_log_posterior, data = fit)

  testthat::expect_type(direct_vals, "double")
  testthat::expect_length(direct_vals, length(take))
  testthat::expect_equal(
    unname(direct_vals),
    unname(lp_vec_all[take]),
    tolerance = 1e-4
  )

  # Basic input validation
  testthat::expect_error(.cmdstan_log_posterior(fit = "not-a-fit", data = list()),
                         info = "Helper should reject invalid input when mis-called with named args.")
  # Proper signature rejects badly shaped x
  testthat::expect_error(.cmdstan_log_posterior(x = numeric(), data = fit))
})
