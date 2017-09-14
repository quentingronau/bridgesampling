
context('basic bridge sampling behavior normal')

test_that("bridge sampler matches anlytical value normal example", {

  # library(bridgesampling)
  library(mvtnorm)

  x <- rmvnorm(1e4, mean = rep(0, 2), sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data) {
    -.5*t(s)%*%s
  }
  assign(x = "log_density", value = log_density, envir = .GlobalEnv)

  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)

  # check repetitions > 1
  bridge_normal <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE, repetitions = 2)
  bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_density,
                                 data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE, repetitions = 2)
  bridge_normal_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                    data = NULL, lb = lb, ub = ub,
                                    method = "normal", silent = TRUE, repetitions = 2)
  bridge_warp3_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                   data = NULL, lb = lb, ub = ub,
                                   method = "warp3", silent = TRUE, repetitions = 2)


  expect_equal(bridge_normal$logml, expected = rep(log(2*pi), length(bridge_normal$logml)), tolerance = 0.01)
  expect_equal(bridge_warp3$logml, expected = rep(log(2*pi), length(bridge_warp3$logml)), tolerance = 0.01)
  expect_equal(bridge_normal_c$logml, expected = rep(log(2*pi), length(bridge_normal_c$logml)), tolerance = 0.01)
  expect_equal(bridge_warp3_c$logml, expected = rep(log(2*pi), length(bridge_warp3_c$logml)), tolerance = 0.01)

  expect_equal(bf(bridge_normal, bridge_warp3)$bf, expected = rep(1, 2), tolerance = 0.1)

  # check repetitions = 1
  bridge_normal <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE)
  bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_density,
                                 data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE)
  bridge_normal_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                    data = NULL, lb = lb, ub = ub,
                                    method = "normal", silent = TRUE)
  bridge_warp3_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                   data = NULL, lb = lb, ub = ub,
                                   method = "warp3", silent = TRUE)

  expect_equal(bridge_normal$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_normal_c$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3_c$logml, expected = log(2*pi), tolerance = 0.01)


  # check using dots repetitions > 1
  mu <- c(1, 2)
  x <- rmvnorm(1e4, mean = mu, sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data, ...) {
    -.5*t(s - mu) %*% (s - mu)
  }
  assign(x = "log_density", value = log_density, envir = .GlobalEnv)
  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)
  bridge_normal_dots <- bridge_sampler(samples = x, log_posterior = log_density, mu,
                                       data = NULL, lb = lb, ub = ub, method = "normal",
                                       silent = TRUE, repetitions = 2)
  bridge_warp3_dots <- bridge_sampler(samples = x, log_posterior = log_density, mu,
                                      data = NULL, lb = lb, ub = ub, method = "normal",
                                      silent = TRUE, repetitions = 2)
  bridge_normal_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                         mu, data = NULL, lb = lb, ub = ub,
                                         method = "normal", silent = TRUE, repetitions = 2)
  bridge_warp3_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                        mu, data = NULL, lb = lb, ub = ub,
                                        method = "warp3", silent = TRUE, repetitions = 2)

  expect_equal(bridge_normal_dots$logml, expected = rep(log(2*pi), length(bridge_normal_dots$logml)), tolerance = 0.01)
  expect_equal(bridge_warp3_dots$logml, expected = rep(log(2*pi), length(bridge_warp3_dots$logml)), tolerance = 0.01)
  expect_equal(bridge_normal_c_dots$logml, expected = rep(log(2*pi), length(bridge_normal_c_dots$logml)), tolerance = 0.01)
  expect_equal(bridge_warp3_c_dots$logml, expected = rep(log(2*pi), length(bridge_warp3_c_dots$logml)), tolerance = 0.01)

  # check using dots
  mu <- c(1, 2)
  x <- rmvnorm(1e4, mean = mu, sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data, ...) {
    -.5*t(s - mu) %*% (s - mu)
  }
  assign(x = "log_density", value = log_density, envir = .GlobalEnv)
  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)
  bridge_normal_dots <- bridge_sampler(samples = x, log_posterior = log_density, mu,
                                       data = NULL, lb = lb, ub = ub, method = "normal",
                                       silent = TRUE)
  bridge_warp3_dots <- bridge_sampler(samples = x, log_posterior = log_density, mu,
                                      data = NULL, lb = lb, ub = ub, method = "normal",
                                      silent = TRUE)
  bridge_normal_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                         mu, data = NULL, lb = lb, ub = ub,
                                         method = "normal", silent = TRUE)
  bridge_warp3_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                        mu, data = NULL, lb = lb, ub = ub,
                                        method = "warp3", silent = TRUE)

  expect_equal(bridge_normal_dots$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3_dots$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_normal_c_dots$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3_c_dots$logml, expected = log(2*pi), tolerance = 0.01)


  # check error_measures
  err <- error_measures(bridge_normal)
  expect_equal(names(err), c("re2", "cv", "percentage"))
  expect_is(unlist(err), "character")

  expect_error(error_measures(bridge_warp3), "not implemented for warp3")

  ### these are meant to check the bf and post_prob functions and not as a meaningful comparisons
  bf <- bf(bridge_normal, bridge_warp3)
  expect_is(bf$bf, "numeric")

  # without prior_prob
  post1 <- post_prob(bridge_normal, bridge_warp3, bridge_normal_c, bridge_warp3_c)
  expect_equal(sum(post1), 1)

  # with prior_prob
  post2 <- post_prob(bridge_normal, bridge_warp3, bridge_normal_c,
                             bridge_warp3_c, prior_prob = c(0.2, 0.1, 0.25, 0.45))
  expect_equal(sum(post2), 1)

  # with incorrect prior_prob
  expect_error(post_prob(bridge_normal, bridge_warp3, bridge_normal_c,
                                 bridge_warp3_c, prior_prob = c(0.2, 0.1, 0.25, 0.55)),
               "do not sum to one")

})
