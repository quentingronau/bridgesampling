#'Computes log marginal likelihood via bridge sampling.
#'@title Log Marginal Likelihood via Bridge Sampling
#'@name bridge_sampler
#'@param samples an \code{mcmc.list} object, a fitted \code{stanfit} object, a
#'  \code{stanreg} object, an \code{rjags} object, a \code{runjags} object, or a
#'  \code{matrix} with posterior samples (\code{colnames} need to correspond to
#'  parameter names in \code{lb} and \code{ub})  with posterior samples.
#'@param log_posterior function or name of function that takes a parameter
#'  vector and the \code{data} as input and returns the log of the unnormalized
#'  posterior density (i.e., a scalar value). If the function name is passed,
#'  the function should exist in the \code{.GlobalEnv}. For special behavior if
#'  \code{cores > 1} see \code{Details}.
#'@param ... additional arguments passed to \code{log_posterior}. Ignored for
#'  the \code{stanfit} and \code{stanreg} methods.
#'@param data data object which is used in \code{log_posterior}.
#'@param stanfit_model for the \code{stanfit} method, an additional object of
#'  class \code{"stanfit"} with the same model as \code{samples}, which will be
#'  used for evaluating the \code{log_posterior} (i.e., it does not need to
#'  contain any samples). The default is to use \code{samples}. In case
#'  \code{samples} was compiled in a different R session or on another computer
#'  with a different OS or setup, the \code{samples} model usually cannot be
#'  used for evaluation. In this case, one can compile the model on the current
#'  computer with \code{iter = 0} and pass it here (this usually needs to be
#'  done before \code{samples} is loaded).
#'@param lb named vector with lower bounds for parameters.
#'@param ub named vector with upper bounds for parameters.
#'@param repetitions number of repetitions.
#'@param method either \code{"normal"} or \code{"warp3"}.
#'@param cores number of cores used for evaluating \code{log_posterior}. On
#'  unix-like systems (where \code{.Platform$OS.type == "unix"} evaluates to
#'  \code{TRUE}; e.g., Linux and Mac OS) forking via \code{\link{mclapply}} is
#'  used. Hence elements needed for evaluation should be in the
#'  \code{\link{.GlobalEnv}}. For other systems (e.g., Windows)
#'  \code{\link{makeCluster}} is used and further arguments specified below will
#'  be used.
#'@param use_neff Boolean which determines whether the effective sample size is
#'  used in the optimal bridge function. Default is TRUE. If FALSE, the number
#'  of samples is used instead. If \code{samples} is a \code{matrix}, it is
#'  assumed that the \code{matrix} contains the samples of one chain in order.
#'  If \code{samples} come from more than one chain, we recommend to use an
#'  \code{mcmc.list} object for optimal performance.
#'@param packages character vector with names of packages needed for evaluating
#'  \code{log_posterior} in parallel (only relevant if \code{cores > 1} and
#'  \code{.Platform$OS.type != "unix"}).
#'@param varlist character vector with names of variables needed for evaluating
#'  \code{log_posterior} (only needed if \code{cores > 1}  and
#'  \code{.Platform$OS.type != "unix"} as these objects will be exported to the
#'  nodes). These objects need to exist in \code{envir}.
#'@param envir specifies the environment for \code{varlist} (only needed if
#'  \code{cores > 1}  and \code{.Platform$OS.type != "unix"} as these objects
#'  will be exported to the nodes). Default is \code{\link{.GlobalEnv}}.
#'@param rcppFile in case \code{cores > 1} and \code{log_posterior} is an
#'  \code{Rcpp} function, \code{rcppFile} specifies the path to the cpp file
#'  (will be compiled on all cores).
#'@param maxiter maximum number of iterations for the iterative updating scheme.
#'  Default is 1,000 to avoid infinite loops.
#'@param param_types character vector of length \code{ncol(samples)} with
#'  \code{"real"}, \code{"simplex"} or \code{"circular"}. For all regular
#'  bounded or unbounded continuous parameters, this should just be
#'  \code{"real"}. However, if there are parameters which lie on a simplex or on
#'  the circle, this should be noted here. Simplex parameters are parameters
#'  which are bounded below by zero and collectively sum to one, such as weights
#'  in a mixture model. For these, the stick-breaking transformation is
#'  performed as described in the Stan reference manual. The circular variables
#'  are given a numerical representation to which the normal distribution is
#'  most likely a good fit. Only possible to use with
#'  \code{bridge_sampler.matrix}.
#'@param silent Boolean which determines whether to print the number of
#'  iterations of the updating scheme to the console. Default is FALSE.
#'@param verbose Boolean. Should internal debug information be printed to
#'  console? Default is \code{FALSE}.
#'@details Bridge sampling is implemented as described in Meng and Wong (1996,
#'  see equation 4.1) using the "optimal" bridge function. When \code{method =
#'  "normal"}, the proposal distribution is a multivariate normal distribution
#'  with mean vector equal to the sample mean vector of \code{samples} and
#'  covariance matrix equal to the sample covariance matrix of \code{samples}.
#'  For a recent tutorial on bridge sampling, see Gronau et al. (in press).
#'
#'  When \code{method = "warp3"}, the proposal distribution is a standard
#'  multivariate normal distribution and the posterior distribution is "warped"
#'  (Meng & Schilling, 2002) so that it has the same mean vector, covariance
#'  matrix, and skew as the samples. \code{method = "warp3"} takes approximately
#'  twice as long as \code{method = "normal"}.
#'
#'  Note that for the \code{matrix} method, the lower and upper bound of a
#'  parameter cannot be a function of the bounds of another parameter.
#'  Furthermore, constraints that depend on multiple parameters of the model are
#'  not supported. This usually excludes, for example, parameters that
#'  constitute a covariance matrix or sets of parameters that need to sum to
#'  one.
#'
#'  However, if the retransformations are part of the model itself and the
#'  \code{log_posterior} accepts parameters on the real line and performs the
#'  appropriate Jacobian adjustments, such as done for \code{stanfit} and
#'  \code{stanreg} objects, such constraints are obviously possible (i.e., we
#'  currently do not know of any parameter supported within Stan that does not
#'  work with the current implementation through a \code{stanfit} object).
#'
#'  \subsection{Parallel Computation}{ On unix-like systems forking is used via
#'  \code{\link{mclapply}}. Hence elements needed for evaluation of
#'  \code{log_posterior} should be in the \code{\link{.GlobalEnv}}.
#'
#'  On other OSes (e.g., Windows), things can get more complicated. For normal
#'  parallel computation, the \code{log_posterior} function can be passed as
#'  both function and function name. If the latter, it needs to exist in the
#'  environment specified in the \code{envir} argument. For parallel computation
#'  when using an \code{Rcpp} function, \code{log_posterior} can only be passed
#'  as the function name (i.e., character). This function needs to result from
#'  calling \code{sourceCpp} on the file specified in \code{rcppFile}.
#'
#'  Due to the way \code{rstan} currently works, parallel computations with
#'  \code{stanfit} and \code{stanreg} objects only work with forking (i.e., NOT
#'  on Windows). }
#'@return if \code{repetitions = 1}, returns a list of class \code{"bridge"}
#'  with components: \itemize{ \item \code{logml}: estimate of log marginal
#'  likelihood. \item \code{niter}: number of iterations of the iterative
#'  updating scheme. \item \code{method}: bridge sampling method that was used
#'  to obtain the estimate. \item \code{q11}: log posterior evaluations for
#'  posterior samples. \item \code{q12}: log proposal evaluations for posterior
#'  samples. \item \code{q21}: log posterior evaluations for samples from
#'  proposal. \item \code{q22}: log proposal evaluations for samples from
#'  proposal. } if \code{repetitions > 1}, returns a list of class
#'  \code{"bridge_list"} with components: \itemize{ \item \code{logml}: numeric
#'  vector with estimates of log marginal likelihood. \item \code{niter}:
#'  numeric vector with number of iterations of the iterative updating scheme
#'  for each repetition. \item \code{method}: bridge sampling method that was
#'  used to obtain the estimates. \item \code{repetitions}: number of
#'  repetitions. }
#'@section Warning: Note that the results depend strongly on the parameter
#'  priors. Therefore, it is strongly advised to think carefully about the
#'  priors before calculating marginal likelihoods. For example, the prior
#'  choices implemented in \pkg{rstanarm} or \pkg{brms} might not be optimal
#'  from a testing point of view. We recommend to use priors that have been
#'  chosen from a testing and not a purely estimation perspective.
#'
#'  Also note that for testing, the number of posterior samples usually needs to
#'  be substantially larger than for estimation.
#'@note To be able to use a \code{stanreg} object for \code{samples}, the user
#'  crucially needs to have specified the \code{diagnostic_file} when fitting
#'  the model in \pkg{rstanarm}.
#'@author Quentin F. Gronau and Henrik Singmann. Parallel computing (i.e.,
#'  \code{cores > 1}) and the \code{stanfit} method use code from \code{rstan}
#'  by Jiaqing Guo, Jonah Gabry, and Ben Goodrich. Ben Goodrich added the
#'  \code{stanreg} method. Kees Mulder added methods for simplex and circular
#'  variables.
#'@references Gronau, Q. F., Sarafoglou, A., Matzke, D., Ly, A., Boehm, U.,
#'  Marsman, M., Leslie, D. S., Forster, J. J., Wagenmakers, E.-J., &
#'  Steingroever, H. (in press). A tutorial on bridge sampling. \emph{Journal of
#'  Mathematical Psychology}. \url{https://arxiv.org/abs/1703.05984} \cr
#'  \code{vignette("bridgesampling_tutorial")}
#'
#'  Gronau, Q. F., Wagenmakers, E.-J., Heck, D. W., & Matzke, D. (2017). \emph{A
#'  simple method for comparing complex models: Bayesian model comparison for
#'  hierarchical multinomial processing tree models using Warp-III bridge
#'  sampling}. Manuscript submitted for publication.
#'  \url{https://psyarxiv.com/yxhfm}
#'
#'  Meng, X.-L., & Wong, W. H. (1996). Simulating ratios of normalizing
#'  constants via a simple identity: A theoretical exploration. \emph{Statistica
#'  Sinica, 6}, 831-860.
#'  \url{http://www3.stat.sinica.edu.tw/statistica/j6n4/j6n43/j6n43.htm}
#'
#'  Meng, X.-L., & Schilling, S. (2002). Warp bridge sampling. \emph{Journal of
#'  Computational and Graphical Statistics, 11(3)}, 552-586.
#'  \url{http://dx.doi.org/10.1198/106186002457}
#'
#'  Overstall, A. M., & Forster, J. J. (2010). Default Bayesian model
#'  determination methods for generalised linear mixed models.
#'  \emph{Computational Statistics & Data Analysis, 54}, 3269-3288.
#'  \url{http://dx.doi.org/10.1016/j.csda.2010.03.008}
#'@example examples/example.bridge_sampler.R
#'
#'@seealso \code{\link{bf}} allows the user to calculate Bayes factors and
#'  \code{\link{post_prob}} allows the user to calculate posterior model
#'  probabilities from bridge sampling estimates. \code{\link{bridge-methods}}
#'  lists some additional methods that automatically invoke the
#'  \code{\link{error_measures}} function.
#'
#'@importFrom mvtnorm rmvnorm dmvnorm
#'@importFrom Matrix nearPD
#'@import Brobdingnag
#'@importFrom stringr str_sub
#'@importFrom stats qnorm pnorm dnorm median cov var
#'@export
bridge_sampler <- function(samples, ...) {
   UseMethod("bridge_sampler", samples)
}

#' @rdname bridge_sampler
#' @export
bridge_sampler.stanfit <- function(samples = NULL, stanfit_model = samples,
                                   repetitions = 1, method = "normal", cores = 1,
                                   use_neff = TRUE, maxiter = 1000, silent = FALSE,
                                   verbose = FALSE, ...) {

  # convert samples into matrix
  if (!requireNamespace("rstan")) stop("package rstan required")
  ex <- rstan::extract(samples, permuted = FALSE)
  skeleton <- .create_skeleton(samples@sim$pars_oi,
                               samples@par_dims[samples@sim$pars_oi])
  upars <- apply(ex, 1:2, FUN = function(theta) {
    rstan::unconstrain_pars(stanfit_model, .rstan_relist(theta, skeleton))
  })

  if (length(dim(upars)) == 2) { # for one parameter models
    dim(upars) <- c(1, dim(upars))
  }

  nr <- dim(upars)[2]
  samples4fit_index <- seq_len(nr) %in% seq_len(round(nr/2)) # split samples in two parts
  samples_4_fit <- apply(upars[,samples4fit_index,,drop=FALSE], 1, rbind)

  samples_4_iter_stan <- upars[,!samples4fit_index,,drop=FALSE]
  samples_4_iter_tmp <- vector("list", dim(upars)[3])
  for (i in seq_along(samples_4_iter_tmp)) {
    samples_4_iter_tmp[[i]] <- coda::as.mcmc(t(samples_4_iter_stan[,,i]))
  }
  samples_4_iter_tmp <- coda::as.mcmc.list(samples_4_iter_tmp)

  if (use_neff) {
    neff <- tryCatch(median(coda::effectiveSize(samples_4_iter_tmp)), error = function(e) {
      warning("effective sample size cannot be calculated, has been replaced by number of samples.", call. = FALSE)
      return(NULL)
    })
  } else {
    neff <- NULL
  }

  samples_4_iter <- apply(samples_4_iter_stan, 1, rbind)

  parameters <- paste0("x", (seq_len(dim(upars)[1])))

  transTypes <- rep("unbounded", length(parameters))
  names(transTypes) <- parameters

  # prepare lb and ub
  lb <- rep(-Inf, length(parameters))
  ub <- rep(Inf, length(parameters))
  names(lb) <- names(ub) <- parameters

  colnames(samples_4_iter) <- paste0("trans_", parameters)
  colnames(samples_4_fit) <- paste0("trans_", parameters)

  # cores > 1 only for unix:
  if (!(.Platform$OS.type == "unix") & (cores != 1)) {
    warning("cores > 1 only possible on Unix/MacOs. Uses 'core = 1' instead.", call. = FALSE)
    cores <- 1L
  }

  # run bridge sampling
  if (cores == 1) {
    bridge_output <- do.call(what = paste0(".bridge.sampler.", method),
                             args = list(samples_4_fit = samples_4_fit,
                                         samples_4_iter = samples_4_iter,
                                         neff = neff,
                                         log_posterior = .stan_log_posterior,
                                         data = list(stanfit = stanfit_model),
                                         lb = lb, ub = ub,
                                         param_types = rep("real", ncol(samples_4_fit)),
                                         transTypes = transTypes,
                                         repetitions = repetitions, cores = cores,
                                         packages = "rstan", maxiter = maxiter, silent = silent,
                                         verbose = verbose,
                                         r0 = 0.5, tol1 = 1e-10, tol2 = 1e-4))
  } else {
    bridge_output <- do.call(what = paste0(".bridge.sampler.", method),
                             args = list(samples_4_fit = samples_4_fit,
                                         samples_4_iter = samples_4_iter,
                                         neff = neff,
                                         log_posterior = .stan_log_posterior,
                                         data = list(stanfit = stanfit_model),
                                         lb = lb, ub = ub,
                                         param_types = rep("real", ncol(samples_4_fit)),
                                         transTypes = transTypes,
                                         repetitions = repetitions, varlist = "stanfit",
                                         envir = sys.frame(sys.nframe()),
                                         cores = cores, packages = "rstan", maxiter = maxiter,
                                         silent = silent, verbose = verbose,
                                         r0 = 0.5, tol1 = 1e-10, tol2 = 1e-4))
  }

  return(bridge_output)

}

#' @rdname bridge_sampler
#' @export
bridge_sampler.mcmc.list <- function(samples = NULL, log_posterior = NULL, ..., data = NULL,
                                     lb = NULL, ub = NULL, repetitions = 1,
                                     param_types = rep("real", ncol(samples[[1]])),
                                     method = "normal", cores = 1, use_neff = TRUE,
                                     packages = NULL, varlist = NULL, envir = .GlobalEnv,
                                     rcppFile = NULL, maxiter = 1000, silent = FALSE,
                                     verbose = FALSE) {
  # split samples in two parts
  nr <- nrow(samples[[1]])
  samples4fit_index <- seq_len(nr) %in% seq_len(round(nr/2))
  samples_4_fit_tmp <- samples[samples4fit_index,,drop=FALSE]
  samples_4_fit_tmp <- do.call("rbind", samples_4_fit_tmp)

  # check lb and ub
  if (!is.numeric(lb))
    stop("lb needs to be numeric", call. = FALSE)
  if (!is.numeric(ub))
    stop("ub needs to be numeric", call. = FALSE)
  if (!all(colnames(samples_4_fit_tmp) %in% names(lb)))
    stop("lb does not contain all parameters.", call. = FALSE)
  if (!all(colnames(samples_4_fit_tmp) %in% names(ub)))
    stop("ub does not contain all parameters.", call. = FALSE)

  # transform parameters to real line
  tmp <- .transform2Real(samples_4_fit_tmp, lb, ub)
  samples_4_fit <- tmp$theta_t
  transTypes <- tmp$transTypes
  samples_4_iter_tmp <- lapply(samples[!samples4fit_index,,drop=FALSE],
                               function(x) .transform2Real(x, lb = lb, ub = ub)$theta_t)

  # compute effective sample size
  if (use_neff) {
    samples_4_iter_tmp <- coda::mcmc.list(lapply(samples_4_iter_tmp, coda::mcmc))
    neff <- tryCatch(median(coda::effectiveSize(samples_4_iter_tmp)), error = function(e) {
      warning("effective sample size cannot be calculated, has been replaced by number of samples.", call. = FALSE)
      return(NULL)
    })
  } else {
    neff <- NULL
  }

  # convert to matrix
  samples_4_iter <- do.call("rbind", samples_4_iter_tmp)

  # run bridge sampling
  out <- do.call(what = paste0(".bridge.sampler.", method),
                 args = list(samples_4_fit = samples_4_fit,
                             samples_4_iter = samples_4_iter,
                             neff = neff,
                             log_posterior = log_posterior,
                             "..." = ..., data = data,
                             lb = lb, ub = ub,
                             transTypes = transTypes,
                             repetitions = repetitions, cores = cores,
                             packages = packages, varlist = varlist, envir = envir,
                             param_types = param_types,
                             rcppFile = rcppFile, maxiter = maxiter,
                             silent = silent, verbose = verbose,
                             r0 = 0.5, tol1 = 1e-10, tol2 = 1e-4))

  return(out)

}

#' @rdname bridge_sampler
#' @export
bridge_sampler.mcmc <- function(samples = NULL, log_posterior = NULL, ...,
                                data = NULL, lb = NULL, ub = NULL,
                                repetitions = 1, method = "normal",
                                cores = 1, use_neff = TRUE,
                                packages = NULL, varlist = NULL,
                                envir = .GlobalEnv, rcppFile = NULL,
                                maxiter = 1000,
                                param_types = rep("real", ncol(samples)),
                                silent = FALSE, verbose = FALSE) {
  samples <- as.matrix(samples)
  bridge_output <- bridge_sampler(samples = samples,
                                  log_posterior = log_posterior,
                                  ...,
                                  data = data, lb = lb, ub = ub,
                                  repetitions = repetitions,
                                  method = method,
                                  cores = cores, use_neff = use_neff,
                                  packages = packages, varlist = varlist,
                                  envir = envir, rcppFile = rcppFile,
                                  maxiter = maxiter,
                                  param_types = param_types,
                                  silent = silent, verbose = verbose)
  return(bridge_output)
}

#' @export
#' @rdname bridge_sampler
bridge_sampler.matrix <- function(samples = NULL, log_posterior = NULL, ...,
                                data = NULL, lb = NULL, ub = NULL,
                                repetitions = 1, method = "normal",
                                cores = 1, use_neff = TRUE,
                                packages = NULL, varlist = NULL,
                                envir = .GlobalEnv, rcppFile = NULL,
                                maxiter = 1000,
                                param_types = rep("real", ncol(samples)),
                                silent = FALSE, verbose = FALSE) {

  # see Meng & Wong (1996), equation 4.1

  # Check simplex computation
  is_simplex_param <- param_types == "simplex"
  if (any(is_simplex_param)) {
    simplex_samples <- samples[, is_simplex_param]

    if (any(!(round(rowSums(simplex_samples), 6) == 1L))) {
      stop(paste("Simplex parameters do not sum to one. This could be due to
               having multiple separate sets of simplex parameters, which are
               not supported. "))
    }

    # Remove the last simplex variable because it is superfluous.
    last_sim <- which(is_simplex_param)[sum(is_simplex_param)]
    samples <- samples[, -last_sim]
    param_types <- param_types[-last_sim]
    lb <- lb[-last_sim]
    ub <- ub[-last_sim]
  }

  # transform parameters to real line
  tmp <- .transform2Real(samples, lb, ub, theta_types = param_types)
  theta_t <- tmp$theta_t
  transTypes <- tmp$transTypes

  # split samples for proposal/iterative scheme
  nr <- nrow(samples)
  samples4fit_index <- seq_len(nr) %in% seq_len(round(nr/2)) # split samples in two parts
  samples_4_fit <- theta_t[samples4fit_index, ,drop = FALSE]
  samples_4_iter <- theta_t[!samples4fit_index, , drop = FALSE]

  # compute effective sample size
  if (use_neff) {
    neff <- tryCatch(median(coda::effectiveSize(coda::mcmc(samples_4_iter))),
                     error = function(e) {
                       warning("effective sample size cannot be calculated, has been replaced by number of samples.", call. = FALSE)
                       return(NULL)
                     })
  } else {
    neff <- NULL
  }

  out <- do.call(what = paste0(".bridge.sampler.", method),
                 args = list(samples_4_fit = samples_4_fit,
                             samples_4_iter = samples_4_iter,
                             neff = neff,
                             log_posterior = log_posterior,
                             "..." = ..., data = data,
                             lb = lb, ub = ub,
                             transTypes = transTypes,
                             param_types = param_types,
                             repetitions = repetitions, cores = cores,
                             packages = packages, varlist = varlist, envir = envir,
                             rcppFile = rcppFile, maxiter = maxiter,
                             silent = silent, verbose = verbose,
                             r0 = 0.5, tol1 = 1e-10, tol2 = 1e-4))
  return(out)

}

#' @rdname bridge_sampler
#' @export
#' @importFrom utils read.csv
bridge_sampler.stanreg <-
  function(samples, repetitions = 1, method = "normal", cores = 1,
           use_neff = TRUE, maxiter = 1000, silent = FALSE,
           verbose = FALSE, ...) {
    df <- eval(samples$call$diagnostic_file)
    if (is.null(df))
      stop("the 'diagnostic_file' option must be specified in the call to ",
           samples$stan_function, " to use the 'bridge_sampler'")
    sf <- samples$stanfit
    chains <- ncol(sf)
    if (chains > 1) df <- sapply(1:chains, FUN = function(j)
      sub("\\.csv$", paste0("_", j, ".csv"), df))
    samples_list <- lapply(df, FUN = function(f) {
      d <- read.csv(f, comment.char = "#")
      excl <- c("lp__", "accept_stat__", "stepsize__" ,"treedepth__",
                "n_leapfrog__", "divergent__", "energy__")
      d <- d[,!(colnames(d) %in% excl), drop = FALSE]
      coda::as.mcmc(as.matrix(d[, 1:rstan::get_num_upars(sf), drop = FALSE]))
    })
    samples <- coda::as.mcmc.list(samples_list)
    lb <- rep(-Inf, ncol(samples[[1]]))
    ub <- rep( Inf, ncol(samples[[1]]))
    names(lb) <- names(ub) <- colnames(samples[[1]])

    # cores > 1 only for unix:
    if (!(.Platform$OS.type == "unix") & (cores != 1)) {
      warning("cores > 1 only possible on Unix/MacOs. Uses 'core = 1' instead.", call. = FALSE)
      cores <- 1L
    }

    if (cores == 1) {
      bridge_output <- bridge_sampler(samples = samples, log_posterior = .stan_log_posterior,
                                      data = list(stanfit = sf), lb = lb, ub = ub,
                                      repetitions = repetitions, method = method, cores = cores,
                                      use_neff = use_neff, packages = "rstan",
                                      maxiter = maxiter, silent = silent,
                                      verbose = verbose)
    } else {
      bridge_output <- bridge_sampler(samples = samples,
                                      log_posterior = .stan_log_posterior,
                                      data = list(stanfit = sf), lb = lb, ub = ub,
                                      repetitions = repetitions, varlist = "stanfit",
                                      envir = sys.frame(sys.nframe()), method = method,
                                      cores = cores, use_neff = use_neff,
                                      packages = "rstan", maxiter = maxiter,
                                      silent = silent, verbose = verbose)
    }
    return(bridge_output)
}

#' @rdname bridge_sampler
#' @export
bridge_sampler.rjags <- function(samples = NULL, log_posterior = NULL, ..., data = NULL,
                                 lb = NULL, ub = NULL, repetitions = 1,
                                 method = "normal", cores = 1, use_neff = TRUE,
                                 packages = NULL, varlist = NULL,
                                 envir = .GlobalEnv, rcppFile = NULL,
                                 maxiter = 1000, silent = FALSE, verbose = FALSE) {


  # convert to mcmc.list
  cn <- colnames(samples$BUGSoutput$sims.matrix)
  samples <- coda::as.mcmc(samples)
  samples <- samples[,cn != "deviance", drop = FALSE]

  # run bridge sampling
  out <- bridge_sampler(samples = samples, log_posterior = log_posterior, ...,
                        data = data, lb = lb, ub = ub, repetitions = repetitions,
                        method = method, cores = cores, use_neff = use_neff,
                        packages = packages, varlist = varlist, envir = envir,
                        rcppFile = rcppFile, maxiter = maxiter, silent = silent,
                        verbose = verbose)

  return(out)

}

#' @rdname bridge_sampler
#' @export
bridge_sampler.runjags <- function(samples = NULL, log_posterior = NULL, ..., data = NULL,
                                   lb = NULL, ub = NULL, repetitions = 1,
                                   method = "normal", cores = 1, use_neff = TRUE,
                                   packages = NULL, varlist = NULL,
                                   envir = .GlobalEnv, rcppFile = NULL,
                                   maxiter = 1000, silent = FALSE, verbose = FALSE) {


  # convert to mcmc.list
  samples <- coda::as.mcmc.list(samples)

  # run bridge sampling
  out <- bridge_sampler(samples = samples, log_posterior = log_posterior, ...,
                        data = data, lb = lb, ub = ub, repetitions = repetitions,
                        method = method, cores = cores, use_neff = use_neff,
                        packages = packages, varlist = varlist, envir = envir,
                        rcppFile = rcppFile, maxiter = maxiter, silent = silent,
                        verbose = verbose)

  return(out)

}

#' @rdname bridge_sampler
#' @export
bridge_sampler.MCMC_refClass <- function(samples,
                                  repetitions = 1,
                                  method = "normal",
                                  cores = 1,
                                  use_neff = TRUE,
                                  maxiter = 1000,
                                  silent = FALSE,
                                  verbose = FALSE,
                                  ...) {
  if (!requireNamespace("nimble")) stop("package nimble required")

  ## functions for nimble support
  .log_posterior_nimble <- ".log_posterior_nimble <- nimble::nimbleFunction(

    # based on code by Perry de Valpine

    ## setup code is executed in R and specializes an instance
    ## of the nimbleFunction to a particular model or nodes
    setup = function(model, nodes) {
      calcNodes <- model$getDependencies(nodes)
    },
    ## run code is called repeatedly and can be converted into C++
    run = function(sample = double(1)) {
      values(model, nodes) <<- sample
      out <- model$calculate(calcNodes)
      return(out)
      returnType(double(0))
    }
  )"
  eval(parse(text = .log_posterior_nimble)) ## trick for avoiding R CMD check NOTEs
  .nimble_bounds <- function(samples, model, which) {

    if ( ! (which %in% c("lower", "upper")) ) {
      stop('"which" needs to be either "lower" or "upper"\n',  call. = FALSE)
    }

    cn <- colnames(samples)
    bounds <- numeric(length(cn))
    names(bounds) <- cn

    for (i in seq_along(cn)) {
      bounds[[cn[i]]] <- model$getBound(cn[i], which)
    }

    return(bounds)

  }

  # cores > 1 only for unix:
  if (!(.Platform$OS.type == "unix") & (cores != 1)) {
    warning("cores > 1 only possible on Unix/MacOs. Uses 'core = 1' instead.",
            call. = FALSE)
    cores <- 1L
  }

  mcmc_samples <- as.matrix(samples$mvSamples)

  if (all(is.na(mcmc_samples))) {
    stop("nimble object does not contain samples. Call runMCMC() first.",
         call. = FALSE)
  }

  # make sure that samples is a list
  if (is.matrix(mcmc_samples)) {
    # TRUE in case nchains = 1
    mcmc_samples <- list(mcmc_samples)
  }

  # convert samples to mcmc.list
  samples_mcmc <- lapply(mcmc_samples, FUN = coda::as.mcmc)
  samples_mcmc_list <- coda::as.mcmc.list(samples_mcmc)

  ## get model name from MCMC_refClass object
  mod_name <- ls(samples$nimbleProject$models)[1]
  nimble_model <- samples$nimbleProject$models[[mod_name]]

  # compile log_posterior for bridge sampling
  log_posterior_tmp <- .log_posterior_nimble(model = nimble_model,
                                             nodes = colnames(mcmc_samples[[1]]))
  suppressMessages(
    clog_posterior <- nimble::compileNimble(log_posterior_tmp,
                                            project = nimble_model))

  # wrapper to match required format for log_posterior
  log_posterior <- function(x, data) {
    clog_posterior$run(x)
  }

  out <- bridge_sampler(samples = samples_mcmc_list,
                        log_posterior = log_posterior,
                        ...,
                        data = NULL,
                        lb = .nimble_bounds(mcmc_samples[[1]],
                                            nimble_model, "lower"),
                        ub = .nimble_bounds(mcmc_samples[[1]],
                                            nimble_model, "upper"),
                        repetitions = repetitions,
                        method = method,
                        cores = cores,
                        use_neff = use_neff,
                        packages = "nimble",
                        maxiter = maxiter,
                        silent = silent,
                        verbose = verbose)

  return(out)

}

