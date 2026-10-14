#--------------------------------------------------------------------------
# functions for Stan support via rstan
#--------------------------------------------------------------------------

# taken from rstan:
.rstan_relist <- function(x, skeleton) {
  lst <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) {
    dim(lst[[i]]) <- dim(skeleton[[i]])
  }
  lst
}

# taken from rstan:
.create_skeleton <- function(pars, dims) {
  lst <- lapply(seq_along(pars), function(i) {
    len_dims <- length(dims[[i]])
    if (len_dims < 1) {
      return(0)
    }
    return(array(0, dim = dims[[i]]))
  })
  names(lst) <- pars
  lst
}

.stan_log_posterior <- function(s.row, data) {
  out <- tryCatch(
    rstan::log_prob(object = data$stanfit, upars = s.row),
    error = function(e) -Inf
  )
  if (is.na(out)) {
    out <- -Inf
  }
  return(out)
}

.cmdstan_log_posterior <- function(s.row, data) {
  if ("lp__" %in% names(s.row)) {
    s.row <- s.row[!names(s.row) %in% "lp__"]
  }

  if (!is.numeric(s.row)) {
    s.row <- as.numeric(s.row)
  }
  out <- tryCatch(
    {
      log_prob <- data$log_prob(s.row, jacobian = TRUE)
      log_prob
    },
    error = function(e) {
      print(e)
      -Inf
    }
  )

  if (is.na(out)) {
    out <- -Inf
  }
  return(out)
}

#--------------------------------------------------------------------------
# posterior scores for the score-matching proposal fit
#--------------------------------------------------------------------------

# Evaluate the posterior score s(theta) = grad log p(theta) over a whole
# matrix of draws on the unconstrained scale.
#
# Neither rstan::grad_log_prob() nor cmdstanr's $grad_log_prob() accepts
# more than one draw at a time, so the loop below is unavoidable.
#
# `grad_fun` takes a single unconstrained parameter vector and returns a
# numeric vector of the same length. Rows for which it fails, or for
# which it returns a non-finite value, are filled with NA; they are
# dropped, with a warning, by .bridge.sampler.normal(). Returns a matrix
# with the dimensions of `draws`, or NULL if no row could be evaluated,
# in which case the caller falls back to the sample covariance.
.gradients_at_draws <- function(draws, grad_fun, what = "grad_log_prob()") {
  p <- ncol(draws)
  gm <- matrix(
    NA_real_,
    nrow = nrow(draws),
    ncol = p,
    dimnames = list(NULL, colnames(draws))
  )
  first_error <- NULL
  for (i in seq_len(nrow(draws))) {
    # A failure is a property of the draw, so the row is left as NA and
    # the fit continues.
    g <- tryCatch(
      as.numeric(grad_fun(draws[i, ])),
      error = function(e) {
        if (is.null(first_error)) {
          first_error <<- conditionMessage(e)
        }
        NULL
      }
    )
    if (is.null(g)) {
      next
    }
    # A wrong length is structural and would silently corrupt the fit.
    if (length(g) != p) {
      stop(
        sprintf(
          "%s returned %d value(s) at draw %d but the model has %d parameter(s).",
          what,
          length(g),
          i,
          p
        ),
        call. = FALSE
      )
    }
    gm[i, ] <- g
  }

  if (all(is.na(gm))) {
    warning(
      what,
      " could not be evaluated at any draw",
      if (!is.null(first_error)) paste0(" (", first_error, ")") else "",
      "; falling back to the sample covariance.",
      call. = FALSE
    )
    return(NULL)
  }
  gm
}

# Posterior scores for a cmdstanr fit. Returns NULL when `proposal_fit`
# does not need gradients, or when the model methods that expose
# grad_log_prob() are unavailable.
.cmdstan_gradients <- function(samples, draws, proposal_fit) {
  if (proposal_fit == "sample") {
    return(NULL)
  }
  # grad_log_prob() is exposed by cmdstanr's init_model_methods().
  ok <- tryCatch(
    {
      suppressMessages(samples$init_model_methods())
      TRUE
    },
    error = function(e) FALSE
  )
  if (!ok) {
    warning(
      "could not initialise the cmdstanr model methods needed for ",
      "proposal_fit = 'hybrid'; falling back to the sample covariance.",
      call. = FALSE
    )
    return(NULL)
  }
  .gradients_at_draws(
    draws,
    function(upars) samples$grad_log_prob(unconstrained_variables = upars),
    what = "grad_log_prob()"
  )
}

# rstan counterpart of .cmdstan_gradients(). `upars` holds draws on the
# unconstrained scale, so adjust_transform = TRUE returns the score of
# the same density that .stan_log_posterior() evaluates.
.rstan_gradients <- function(stanfit_model, upars, proposal_fit) {
  if (proposal_fit == "sample") {
    return(NULL)
  }
  .gradients_at_draws(
    upars,
    function(u) {
      rstan::grad_log_prob(
        object = stanfit_model,
        upars = u,
        adjust_transform = TRUE
      )
    },
    what = "rstan::grad_log_prob()"
  )
}
