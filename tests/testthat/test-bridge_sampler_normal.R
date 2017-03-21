
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

})
