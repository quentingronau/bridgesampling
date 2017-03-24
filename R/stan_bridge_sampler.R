
# taken fron rstan
.rstan_relist <- function (x, skeleton) {
  lst <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) dim(lst[[i]]) <- dim(skeleton[[i]])
  lst
}

# taken from rstan
.create_skeleton <- function (pars, dims) {
  lst <- lapply(seq_along(pars), function(i) {
    len_dims <- length(dims[[i]])
    if (len_dims < 1)
      return(0)
    return(array(0, dim = dims[[i]]))
  })
  names(lst) <- pars
  lst
}

.stan_log_posterior_local <- function(s.row, data) {
  rstan::log_prob(object = data$stanfit, upars = s.row)
}

.stan_log_posterior_distributed <- function(s.row, data) {
  rstan::log_prob(object = .stanfit_object_bridgesampling, upars = s.row)
}

.create_stan_log_posterior <- function(stanfit) {
  assign(".stanfit_object_bridgesampling", rstan::sampling(stanfit@stanmodel), envir = .GlobalEnv)
}

#' Computes marginal likelihood for an object of class \code{stanfit}.
#' @export
#' @title Marginal likelihood for an object of class \code{stanfit}
#' @param stanfit an object of class \code{"stanfit"}.
#' @param method either \code{"normal"} or \code{"warp3"}.
#' @param cores number of cores used for computations.
#' @param maxiter maximum number of iterations for the iterative updating scheme. Default is 1,000 to avoid infinite loops.
#' @param silent Boolean which determines whether to print the number of iterations of the updating scheme to the console. Default is FALSE.
#' @details
#' @return a list of class \code{"bridge"} with components:
#' \itemize{
#'  \item \code{logml}: estimate of log marginal likelihood.
#'  \item \code{niter}: number of iterations of the iterative updating scheme.
#'  \item \code{method}: bridge sampling method that was used to obtain the estimate.
#'  \item \code{q11}: log_posterior evaluations for posterior samples.
#'  \item \code{q12}: log proposal evaluations for posterior samples.
#'  \item \code{q21}: log_posterior evaluations for samples from proposal.
#'  \item \code{q22}: log proposal evaluations for samples from proposal.
#' }
#' @author Quentin F. Gronau
#' @import rstan
stan_bridge_sampler <- function(stanfit = NULL, method = "normal", cores = 1,
                                maxiter = 1000, silent = FALSE) {

  # convert samples into matrix
  ex <- extract(stanfit, permuted = FALSE)
  skeleton <- .create_skeleton(stanfit@model_pars, stanfit@par_dims)
  upars <- apply(ex, 1:2, FUN = function(theta) {
    rstan::unconstrain_pars(stanfit, .rstan_relist(theta, skeleton))
  })
  samples <- apply(upars, 1, rbind)

  # prepare lb and ub
  colnames(samples) <- paste0("x", seq_len(ncol(samples)))
  lb <- rep(-Inf, ncol(samples))
  ub <- rep(Inf, ncol(samples))
  names(lb) <- names(ub) <- colnames(samples)

  #browser()
  # run bridge sampling
  if (cores == 1) {
    bridge_output <- bridge_sampler(samples = samples, log_posterior = .stan_log_posterior,
                                    data = list(stanfit = stanfit), lb = lb, ub = ub,
                                    method = method, cores = cores, packages = "rstan",
                                      maxiter = maxiter, silent = silent)
  } else {
    bridge_output <- bridge_sampler(samples = samples,
                                    log_posterior = .stan_log_posterior_distributed,
                                    data = NULL, lb = lb, ub = ub,
                                    varlist = "stanfit", envir = sys.frame(sys.nframe()),
                                    method = method, cores = cores, packages = "rstan",
                                    eval_q =  quote(assign(".stanfit_object_bridgesampling",
                                                     rstan::sampling(stanfit@stanmodel, chains = 0),
                                                     envir = .GlobalEnv)),
                                    maxiter = maxiter, silent = silent)
  }

  return(bridge_output)

}
