
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

.stan_log_posterior <- function(s.row, data) {
  rstan::log_prob(object = data$stanfit, upars = s.row)
}


# .stan_log_posterior_local <- function(s.row, data) {
#   rstan::log_prob(object = data$stanfit, upars = s.row)
# }
#
# .stan_log_posterior_distributed <- function(s.row, data) {
#   rstan::log_prob(object = .stanfit_object_bridgesampling, upars = s.row)
# }
#
# .create_stan_log_posterior <- function(stanfit) {
#   assign(".stanfit_object_bridgesampling", rstan::sampling(stanfit@stanmodel), envir = .GlobalEnv)
#
# }

#' Computes marginal likelihood for an object of class \code{stanfit}.
#' @export
#' @title Marginal likelihood for an object of class \code{stanfit}
#' @param stanfit_data an object of class \code{"stanfit"} that contains the posterior samples.
#' @param stanfit_model an object of class \code{"stanfit"} with the same model as \code{stanfit_data}, which will be used for evaluating the \code{log_posterior} (i.e., it does not need to contain any samples). The default is to use \code{stanfit_data}. In case \code{stanfit_data} was compiled on another computer with a different OS or setup, \code{stanfit_data} usually cannot be used for evaluation. In this case, one can compile the model on the current computer with \code{iter = 0} and pass it here (this usually needs to be done before \code{stanfit_data} is loaded).
#' @param method either \code{"normal"} or \code{"warp3"}.
#' @param cores number of cores used for computations. \code{cores > 1} is currently only supported on unix-like systems that support forking via \code{\link{mclapply}} (e.g., Linux or Mac OS).
#' @param maxiter maximum number of iterations for the iterative updating scheme. Default is 1,000 to avoid infinite loops.
#' @param silent Boolean which determines whether to print the number of iterations of the updating scheme to the console. Default is FALSE.
#' @details Works.
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
#' @author Quentin F. Gronau and Henrik Singmann. Uses code from \code{rstan} by Jiaqing Guo, Jonah Gabry, and Ben Goodrich.
#' @import rstan
stan_bridge_sampler <- function(stanfit_data = NULL, stanfit_model = stanfit_data,
                                method = "normal", cores = 1,
                                maxiter = 1000, silent = FALSE) {

  # convert samples into matrix
  ex <- extract(stanfit_data, permuted = FALSE)
  skeleton <- .create_skeleton(stanfit_model@model_pars, stanfit_model@par_dims)
  upars <- apply(ex, 1:2, FUN = function(theta) {
    rstan::unconstrain_pars(stanfit_model, .rstan_relist(theta, skeleton))
  })

  if (length(dim(upars)) == 3) {
    samples <- apply(upars, 1, rbind)
  } else if (length(dim(upars)) == 2) {
    samples <- as.matrix(as.vector(upars))
  }

  # prepare lb and ub
  colnames(samples) <- paste0("x", seq_len(ncol(samples)))
  lb <- rep(-Inf, ncol(samples))
  ub <- rep(Inf, ncol(samples))
  names(lb) <- names(ub) <- colnames(samples)

  #browser()
  # run bridge sampling
  if (cores == 1) {
    bridge_output <- bridge_sampler(samples = samples, log_posterior = .stan_log_posterior,
                                    data = list(stanfit = stanfit_model), lb = lb, ub = ub,
                                    method = method, cores = cores, packages = "rstan",
                                    maxiter = maxiter, silent = silent)
  } else {
    bridge_output <- bridge_sampler(samples = samples,
                                    log_posterior = .stan_log_posterior,
                                    data = list(stanfit = stanfit_model), lb = lb, ub = ub,
                                    varlist = "stanfit", envir = sys.frame(sys.nframe()),
                                    method = method, cores = cores, packages = "rstan",
                                    maxiter = maxiter, silent = silent)
  }

  return(bridge_output)

}
