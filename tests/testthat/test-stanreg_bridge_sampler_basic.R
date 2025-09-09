
context('bridge_sampler.stanreg works.')

test_that("stan_bridge_sampler", {

  testthat::skip_on_cran()
  testthat::skip_on_ci()

  if (require(rstanarm)) {

    fit_1 <- stan_glm(mpg ~ wt + qsec + am, data = mtcars,
                     chains = 2, cores = 2, iter = 5000,
                     diagnostic_file = file.path(tempdir(), "df.csv"),
                     refresh = 0)
    bridge_norm <- bridge_sampler(fit_1, silent = TRUE)

    fit_2 <- update(fit_1, formula = . ~ . + cyl, refresh = 0)
    bridge_warp <- bridge_sampler(fit_2, method = "warp3", silent = TRUE)

    expect_true(bridge_norm$logml > bridge_warp$logml)
  }
})


test_that("stan_bridge_sampler in multicore", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  #testthat::skip_on_os("windows")
  if (require(rstanarm)) {
    fit_1 <- stan_glm(mpg ~ wt + qsec + am, data = mtcars,
                     chains = 2, cores = 2, iter = 5000,
                     diagnostic_file = file.path(tempdir(), "df.csv"),
                     refresh = 0)
    bridge_norm <- bridge_sampler(fit_1, cores = 2, silent = TRUE)

    fit_2 <- update(fit_1, formula = . ~ . + cyl, refresh = 0)
    bridge_warp <- bridge_sampler(fit_2, method = "warp3", cores = 2, silent = TRUE)

    expect_true(bridge_norm$logml > bridge_warp$logml)
  }
})
