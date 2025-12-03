test_that(".bs_compute_ess uses posterior when available and enabled", {
  skip_if_not_installed("posterior")

  set.seed(123)
  # Simple correlated chain to make ESS less than n
  draws <- cbind(theta = cumsum(rnorm(2000)))

  fn <- getFromNamespace(".bs_compute_ess", "bridgesampling")

  old_opts <- options(bridgesampling.use_posterior_ess = TRUE)
  on.exit(options(old_opts), add = TRUE)

  neff <- fn(draws, use_neff = TRUE)

  posterior_draws <- posterior::as_draws_matrix(draws)
  expected <- as.numeric(median(posterior::ess_mean(posterior_draws)))

  expect_equal(neff, expected)
})

test_that(".bs_compute_ess falls back to coda when posterior is disabled", {
  set.seed(123)
  draws <- cbind(theta = cumsum(rnorm(2000)))

  fn <- getFromNamespace(".bs_compute_ess", "bridgesampling")

  old_opts <- options(bridgesampling.use_posterior_ess = FALSE)
  on.exit(options(old_opts), add = TRUE)

  neff <- fn(draws, use_neff = TRUE)

  mcmc_obj <- coda::mcmc(draws)
  expected <- as.numeric(median(coda::effectiveSize(mcmc_obj)))

  expect_equal(neff, expected)
})

test_that(".bs_compute_ess returns n when use_neff = FALSE", {
  set.seed(123)
  draws <- cbind(theta = rnorm(100))

  fn <- getFromNamespace(".bs_compute_ess", "bridgesampling")

  neff <- fn(draws, use_neff = FALSE)

  expect_identical(neff, nrow(draws))
})

test_that("bridge_sampler runs with posterior-based ESS when available", {
  skip_if_not_installed("posterior")

  set.seed(123)
  x <- rnorm(100)

  log_post <- function(theta, data) {
    # Standard normal likelihood and prior, up to a constant
    sum(dnorm(data, mean = theta, log = TRUE)) + dnorm(theta, log = TRUE)
  }

  samples <- matrix(rnorm(2000), ncol = 1)
  colnames(samples) <- "theta"

  old_opts <- options(bridgesampling.use_posterior_ess = TRUE)
  on.exit(options(old_opts), add = TRUE)

  fit <- bridge_sampler(
    samples      = samples,
    log_posterior = log_post,
    data         = x,
    lb           = -Inf,
    ub           = Inf,
    use_neff     = TRUE,
    method       = "normal",
    silent       = TRUE
  )

  expect_true(is.finite(logml(fit)))
})
