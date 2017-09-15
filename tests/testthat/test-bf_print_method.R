
context('bridge sampling print method')

test_that("bf print method correctly displayed", {

  # library(bridgesampling)
  library(mvtnorm)

  x <- rmvnorm(1e4, mean = rep(0, 2), sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data) {
    -.5*t(s)%*%s
  }

  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)

  # repetitions = 1
  bridge_normal <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE)
  bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE)
  BF <- bf(bridge_normal, bridge_warp3)
  log_BF <- bf(bridge_normal, bridge_warp3, log = TRUE)

  expect_output(print(BF), "The estimated Bayes factor")
  expect_output(print(log_BF), "The estimated log Bayes factor")

  BF2 <- bayes_factor(bridge_normal, bridge_warp3)
  log_BF2 <- bayes_factor(bridge_normal, bridge_warp3, log = TRUE)

  expect_output(print(BF2), "The estimated Bayes factor")
  expect_output(print(log_BF2), "The estimated log Bayes factor")


  # repetitions > 1
  bridge_normal <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE, repetitions = 2)
  bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_density,
                                 data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE, repetitions = 2)

  BF <- bf(bridge_normal, bridge_warp3)
  log_BF <- bf(bridge_normal, bridge_warp3, log = TRUE)

  expect_output(print(BF), "based on the medians")
  expect_output(print(log_BF), "based on the medians")

  # default
  BF <- bf(1, 2)
  log_BF <- bf(1, 2, log = TRUE)

  expect_output(print(BF), "The Bayes factor")
  expect_output(print(log_BF), "The log Bayes factor")

})
