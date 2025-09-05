context(".cmdstan_log_posterior helper function for CmdStanMCMC")

test_that(".cmdstan_log_posterior and bridge_sampler agree with analytical results", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if_not_installed("bridgesampling")

  if (!file.exists(cmdstanr::cmdstan_path())) {
    skip("CmdStan is not installed in the expected path for cmdstanr.")
  }

  # Access the internal helper from bridgesampling
  pkg <- "bridgesampling"
  expect_true(
    exists(".cmdstan_log_posterior", envir = asNamespace(pkg), inherits = FALSE),
    info = "Internal function .cmdstan_log_posterior not found in 'bridgesampling'."
  )
  .cmdstan_log_posterior <- get(".cmdstan_log_posterior", envir = asNamespace(pkg))


  # Normal - Normal Data Model
  set.seed(321)
  N     <- 40L
  sigma <- 1
  y     <- rnorm(N, 0.25, sigma)
  data_list <- list(N = N, y = y, sigma = sigma)

  stan_code <- "
  data { int<lower=1> N; vector[N] y; real<lower=0> sigma; }
  parameters { real mu; }
  model {
    mu ~ normal(0, 1);
    y ~ normal(mu, sigma);
  }"

  tf <- withr::local_tempfile(fileext = ".stan")
  writeLines(stan_code, tf)

  mod <- cmdstanr::cmdstan_model(tf, quiet = TRUE, force_recompile = TRUE)

  fit <- mod$sample(
    data = data_list,
    seed = 404,
    chains = 20,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 10000,
    refresh = 0
  )

  # 1) lp__ vs analytical log posterior (up to a constant)
  # Prior: mu ~ N(0,1); likelihood: y_i | mu ~ N(mu, sigma^2)
  logpost_mu <- function(mu, y, sigma) {
    -0.5 * mu^2 - sum((y - mu)^2) / (2 * sigma^2)
  }

  # Extract draws
  lp_df <- fit$draws(variables = "lp__", format = "df")
  expect_true(nrow(lp_df) > 10)
  lp_vec_all <- lp_df$lp__

  mu_df <- fit$draws(variables = "mu", format = "df")
  expect_equal(nrow(mu_df), length(lp_vec_all))
  mu_vec_all <- mu_df$mu

  lp_anal_all <- vapply(mu_vec_all, logpost_mu, numeric(1), y = y, sigma = sigma)

  # Compare up to additive constant: center both and compare
  lp_cent <- lp_vec_all - mean(lp_vec_all)
  anal_cent <- lp_anal_all - mean(lp_anal_all)

  expect_equal(lp_cent, anal_cent, tolerance = 1e-4)
  expect_gt(stats::cor(lp_vec_all, lp_anal_all), 0.99) # sanity check on correlation

  # 2) bridgesampling's internal q11 agrees with analytical log posterior
  bs_out <- bridgesampling::bridge_sampler(
    samples     = fit,
    repetitions = 1,
    method      = "normal",
    cores       = 1L,
    use_neff    = FALSE,
    silent      = TRUE,
    verbose     = FALSE
  )

  n_all   <- length(lp_vec_all)
  n_fit   <- round(n_all / 2)
  idx_iter <- seq.int(n_fit + 1L, n_all)

  expect_equal(length(bs_out$q11), length(idx_iter))

  # Compute analytical log posterior for the same subset of draws
  mu_iter       <- mu_vec_all[idx_iter]
  lp_anal_iter  <- vapply(mu_iter, logpost_mu, numeric(1), y = y, sigma = sigma)

  # Compare up to an additive constant by centering both
  q11_cent        <- bs_out$q11 - mean(bs_out$q11)
  anal_iter_cent  <- lp_anal_iter - mean(lp_anal_iter)

  expect_equal(unname(q11_cent), unname(anal_iter_cent), tolerance = 1e-4)
  expect_gt(stats::cor(bs_out$q11, lp_anal_iter), 0.99) # additional sanity check


  # 3) Spot-check the internal helper against lp__ directly
  upars <- fit$unconstrain_draws(format = "matrix")
  expect_equal(nrow(upars), length(lp_vec_all))

  take <- seq_len(min(100L, nrow(upars)))
  direct_vals <- apply(upars[take, , drop = FALSE], 1, .cmdstan_log_posterior, data = fit)

  expect_type(direct_vals, "double")
  expect_length(direct_vals, length(take))
  expect_equal(unname(direct_vals), unname(lp_vec_all[take]), tolerance = 1e-4)

  # 4) Basic input validation on the helper
  expect_error(
    .cmdstan_log_posterior(fit = "not-a-fit", data = list()),
    info = "Helper should reject invalid input when mis-called with named args."
  )

  expect_error(
    .cmdstan_log_posterior(x = numeric(), data = fit),
    info = "Proper signature should reject badly shaped x."
  )
})
