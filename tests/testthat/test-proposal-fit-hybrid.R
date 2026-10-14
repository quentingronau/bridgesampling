context("score-matching proposal fit (proposal_fit = 'hybrid')")

# A zero-mean multivariate normal target with known normalizing constant.
# The unnormalized log density is -0.5 * theta' Sigma^-1 theta, so the log
# marginal likelihood is (p/2) log(2 pi) + 0.5 log det Sigma, and the score
# is -Sigma^-1 theta.
gaussian_target <- function(p = 4, n = 4000, seed = 1) {
  set.seed(seed)
  A <- matrix(rnorm(p * p), p, p)
  Sigma <- crossprod(A) + diag(p)
  Sigma_inv <- solve(Sigma)
  draws <- mvtnorm::rmvnorm(n, sigma = Sigma)
  colnames(draws) <- paste0("x", seq_len(p))
  lb <- rep(-Inf, p)
  ub <- rep(Inf, p)
  names(lb) <- names(ub) <- colnames(draws)
  list(
    draws = draws,
    gradients = -draws %*% Sigma_inv,
    log_posterior = function(theta, data) {
      -0.5 * drop(theta %*% Sigma_inv %*% theta)
    },
    logml = (p / 2) * log(2 * pi) + 0.5 * determinant(Sigma, TRUE)$modulus[1],
    lb = lb,
    ub = ub
  )
}

run_bridge <- function(tg, ...) {
  bridge_sampler(
    samples = tg$draws,
    log_posterior = tg$log_posterior,
    data = NULL,
    lb = tg$lb,
    ub = tg$ub,
    silent = TRUE,
    ...
  )
}

test_that("proposal_fit = 'sample' is the default and is recorded as such", {
  tg <- gaussian_target()
  out <- run_bridge(tg)
  expect_s3_class(out, "bridge")
  expect_equal(out$proposal_fit_info$proposal_fit, "sample")
  expect_equal(out$proposal_fit_info$combiner, "none")
  expect_equal(out$logml, tg$logml, tolerance = 0.01)
})

test_that("proposal_fit = 'hybrid' uses the geometric combiner and is accurate", {
  tg <- gaussian_target()
  out <- run_bridge(tg, gradients = tg$gradients, proposal_fit = "hybrid")
  expect_equal(out$proposal_fit_info$proposal_fit, "hybrid")
  expect_equal(out$proposal_fit_info$combiner, "geometric")
  # half of the draws are used to fit the proposal
  expect_equal(out$proposal_fit_info$n_gradients_used, nrow(tg$draws) / 2)
  expect_equal(out$logml, tg$logml, tolerance = 0.01)
})

test_that("proposal_fit = 'hybrid' without gradients falls back silently", {
  tg <- gaussian_target()
  out <- run_bridge(tg, proposal_fit = "hybrid")
  expect_equal(out$proposal_fit_info$combiner, "none")
  expect_equal(out$logml, tg$logml, tolerance = 0.01)
})

test_that("gradients are ignored for bounded parameters, with a warning", {
  tg <- gaussian_target()
  lb <- tg$lb
  ub <- tg$ub
  lb[1] <- -50
  ub[1] <- 50
  expect_warning(
    out <- bridge_sampler(
      samples = tg$draws,
      log_posterior = tg$log_posterior,
      data = NULL,
      lb = lb,
      ub = ub,
      silent = TRUE,
      gradients = tg$gradients,
      proposal_fit = "hybrid"
    ),
    "bounded parameters"
  )
  expect_equal(out$proposal_fit_info$combiner, "none")
})

test_that("mismatched gradient dimensions warn and fall back", {
  tg <- gaussian_target()
  expect_warning(
    out <- run_bridge(
      tg,
      gradients = tg$gradients[, -1, drop = FALSE],
      proposal_fit = "hybrid"
    ),
    "do not match"
  )
  expect_equal(out$proposal_fit_info$combiner, "none")
})

test_that("non-finite gradients warn but the remaining rows are used", {
  tg <- gaussian_target()
  gradients <- tg$gradients
  gradients[1:3, 1] <- NA_real_
  expect_warning(
    out <- run_bridge(tg, gradients = gradients, proposal_fit = "hybrid"),
    "non-finite gradients"
  )
  expect_equal(out$proposal_fit_info$combiner, "geometric")
  expect_equal(out$proposal_fit_info$n_gradients_used, nrow(tg$draws) / 2 - 3)
})

test_that("an unknown proposal_fit is rejected", {
  tg <- gaussian_target()
  expect_error(run_bridge(tg, proposal_fit = "ledoit-wolf"), "should be one of")
})

test_that("warp3 warns that it ignores proposal_fit", {
  tg <- gaussian_target(p = 2, n = 1000)
  expect_warning(
    out <- run_bridge(
      tg,
      method = "warp3",
      gradients = tg$gradients,
      proposal_fit = "hybrid"
    ),
    "only implemented for method"
  )
  expect_equal(out$method, "warp3")
})
