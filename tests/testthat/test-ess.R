test_that(".bs_compute_ess uses posterior when available and selected", {
  skip_if_not_installed("posterior")
  skip_if_not_installed("coda")

  set.seed(123)
  # Simple correlated chain to make ESS less than n
  draws <- cbind(theta = cumsum(rnorm(2000)))

  fn <- getFromNamespace(".bs_compute_ess", "bridgesampling")

  old_opts <- options(bridgesampling.ess_function = "posterior")
  on.exit(options(old_opts), add = TRUE)

  ess <- fn(draws, use_ess = TRUE)

  posterior_draws <- posterior::as_draws_matrix(draws)
  expected <- as.numeric(median(posterior::ess_mean(posterior_draws)))

  expect_equal(ess, expected)
})

test_that(".bs_compute_ess uses coda when coda is selected (even if posterior is installed)", {
  skip_if_not_installed("coda")

  set.seed(123)
  draws <- cbind(theta = cumsum(rnorm(2000)))

  fn <- getFromNamespace(".bs_compute_ess", "bridgesampling")

  old_opts <- options(bridgesampling.ess_function = "coda")
  on.exit(options(old_opts), add = TRUE)

  ess <- fn(draws, use_ess = TRUE)

  mcmc_obj <- coda::mcmc(draws)
  expected <- as.numeric(median(coda::effectiveSize(mcmc_obj)))

  expect_equal(ess, expected)
})

test_that(".bs_compute_ess falls back to coda when option is invalid", {
  skip_if_not_installed("coda")

  set.seed(123)
  draws <- cbind(theta = cumsum(rnorm(2000)))

  fn <- getFromNamespace(".bs_compute_ess", "bridgesampling")

  old_opts <- options(bridgesampling.ess_function = "not-a-real-method")
  on.exit(options(old_opts), add = TRUE)

  ess <- fn(draws, use_ess = TRUE)

  mcmc_obj <- coda::mcmc(draws)
  expected <- as.numeric(median(coda::effectiveSize(mcmc_obj)))

  expect_equal(ess, expected)
})

test_that(".bs_compute_ess returns n when use_ess = FALSE", {
  set.seed(123)
  draws <- cbind(theta = rnorm(100))

  fn <- getFromNamespace(".bs_compute_ess", "bridgesampling")

  ess <- fn(draws, use_ess = FALSE)

  expect_identical(ess, nrow(draws))
})

test_that("bridge_sampler runs with posterior-based ESS when available", {
  skip_if_not_installed("posterior")
  skip_if_not_installed("coda")

  set.seed(123)
  x <- rnorm(100)

  log_post <- function(theta, data) {
    # Standard normal likelihood and prior, up to a constant
    sum(dnorm(data, mean = theta, log = TRUE)) + dnorm(theta, log = TRUE)
  }

  samples <- matrix(rnorm(2000), ncol = 1)
  colnames(samples) <- "theta"

  old_opts <- options(bridgesampling.ess_function = "posterior")
  on.exit(options(old_opts), add = TRUE)

  fit <- bridge_sampler(
    samples        = samples,
    log_posterior  = log_post,
    data           = x,
    lb             = -Inf,
    ub             = Inf,
    use_ess        = TRUE,
    method         = "normal",
    silent         = TRUE
  )

  expect_true(is.finite(logml(fit)))
})
