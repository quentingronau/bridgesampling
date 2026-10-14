context("posterior scores for the score-matching proposal fit")

draws_3x2 <- function() {
  matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    ncol = 2,
    dimnames = list(NULL, c("trans_a", "trans_b"))
  )
}

test_that("every row is evaluated and the column names are kept", {
  draws <- draws_3x2()
  gm <- .gradients_at_draws(draws, function(u) -u)
  expect_true(is.matrix(gm))
  expect_equal(dim(gm), dim(draws))
  expect_equal(colnames(gm), colnames(draws))
  expect_equal(gm, -draws, check.attributes = FALSE)
})

test_that("the gradient function sees one draw at a time", {
  draws <- draws_3x2()
  seen <- list()
  .gradients_at_draws(draws, function(u) {
    seen[[length(seen) + 1L]] <<- u
    rep(0, length(u))
  })
  expect_length(seen, nrow(draws))
  expect_true(all(vapply(seen, length, 1L) == ncol(draws)))
  expect_equal(seen[[2]], draws[2, ], check.attributes = FALSE)
})

test_that("a draw whose gradient errors becomes an NA row", {
  draws <- draws_3x2()
  gm <- .gradients_at_draws(draws, function(u) {
    if (u[1] == 2) stop("no gradient here")
    -u
  })
  expect_true(all(is.na(gm[2, ])))
  expect_equal(gm[-2, ], -draws[-2, ], check.attributes = FALSE)
})

test_that("non-finite gradients are passed through for the fit to drop", {
  draws <- draws_3x2()
  gm <- .gradients_at_draws(draws, function(u) {
    if (u[1] == 3) c(NaN, Inf) else -u
  })
  expect_false(all(is.finite(gm[3, ])))
  expect_true(all(is.finite(gm[1, ])))
})

test_that("failure at every draw returns NULL with one warning", {
  draws <- draws_3x2()
  expect_warning(
    gm <- .gradients_at_draws(
      draws,
      function(u) stop("boom"),
      what = "grad_log_prob()"
    ),
    "could not be evaluated at any draw"
  )
  expect_null(gm)
  # the underlying message is reported once, to help diagnose the cause
  expect_warning(
    .gradients_at_draws(draws, function(u) stop("boom")),
    "boom"
  )
})

test_that("a gradient of the wrong length is a fatal error", {
  draws <- draws_3x2()
  expect_error(
    .gradients_at_draws(draws, function(u) 1, what = "grad_log_prob()"),
    "returned 1 value\\(s\\) at draw 1 but the model has 2 parameter\\(s\\)"
  )
})

test_that("a single-parameter model works", {
  draws <- matrix(c(-1, 0, 1), ncol = 1, dimnames = list(NULL, "trans_mu"))
  gm <- .gradients_at_draws(draws, function(u) -u)
  expect_equal(dim(gm), c(3L, 1L))
  expect_equal(gm, -draws, check.attributes = FALSE)
})

# --- backend wrappers -------------------------------------------------

test_that("no gradients are computed when proposal_fit is 'sample'", {
  draws <- draws_3x2()
  # a stub that would error if it were ever touched
  stub <- list(
    init_model_methods = function(...) stop("must not be called"),
    grad_log_prob = function(...) stop("must not be called")
  )
  expect_null(.cmdstan_gradients(stub, draws, "sample"))
  expect_null(.rstan_gradients(stub, draws, "sample"))
})

test_that(".cmdstan_gradients() drives grad_log_prob() over the draws", {
  draws <- draws_3x2()
  initialised <- 0L
  stub <- list(
    init_model_methods = function(...) {
      initialised <<- initialised + 1L
      invisible(NULL)
    },
    grad_log_prob = function(unconstrained_variables) -unconstrained_variables
  )
  gm <- .cmdstan_gradients(stub, draws, "hybrid")
  expect_equal(gm, -draws, check.attributes = FALSE)
  # the model methods are initialised once, not once per draw
  expect_equal(initialised, 1L)
})

test_that(".cmdstan_gradients() falls back when model methods are unavailable", {
  draws <- draws_3x2()
  stub <- list(
    init_model_methods = function(...) {
      stop("Model methods cannot be used with a pre-compiled Stan executable")
    },
    grad_log_prob = function(...) stop("must not be reached")
  )
  expect_warning(
    gm <- .cmdstan_gradients(stub, draws, "hybrid"),
    "could not initialise the cmdstanr model methods"
  )
  expect_null(gm)
})

test_that("the proposal fit consumes a gradient matrix with NA rows", {
  set.seed(3)
  p <- 3
  n <- 600
  Sigma <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  Sigma_inv <- solve(Sigma)
  draws <- mvtnorm::rmvnorm(n, sigma = Sigma)
  colnames(draws) <- paste0("x", seq_len(p))
  lb <- rep(-Inf, p)
  ub <- rep(Inf, p)
  names(lb) <- names(ub) <- colnames(draws)

  # exactly the shape .gradients_at_draws() returns when two draws fail
  gradients <- .gradients_at_draws(draws, function(u) {
    if (identical(u, draws[1, ]) || identical(u, draws[2, ])) {
      stop("no gradient here")
    }
    -drop(u %*% Sigma_inv)
  })
  expect_true(all(is.na(gradients[1:2, ])))

  expect_warning(
    out <- bridge_sampler(
      samples = draws,
      log_posterior = function(theta, data) {
        -0.5 * drop(theta %*% Sigma_inv %*% theta)
      },
      data = NULL,
      lb = lb,
      ub = ub,
      silent = TRUE,
      gradients = gradients,
      proposal_fit = "hybrid"
    ),
    "2 of 300 draws produced non-finite gradients"
  )
  expect_equal(out$proposal_fit_info$combiner, "geometric")
  expect_equal(out$proposal_fit_info$n_gradients_used, 298)
})
