# library(mvtnorm)
# library(Matrix)
# library(Brobdingnag)
# library(stringr)
# library(coda)

#' @importFrom  stats qnorm
.transform2Real <- function(theta, lb, ub) {

  ### transform samples to real line

  theta_t <- theta
  transTypes <- character()
  cn <- colnames(theta)

  for (i in seq_len(ncol(theta))) {

    p <- cn[i]

    if (lb[[p]] < ub[[p]] && is.infinite(lb[[p]]) && is.infinite(ub[[p]])) {
      transTypes[[p]] <- "unbounded"
      theta_t[,i] <- theta[,i]
    } else if (lb[[p]] < ub[[p]] && is.finite(lb[[p]]) && is.infinite(ub[[p]])) {
      transTypes[[p]] <- "lower-bounded"
      theta_t[,i] <- log(theta[,i] - lb[[p]])
    } else if (lb[[p]] < ub[[p]] && is.infinite(lb[[p]]) && is.finite(ub[[p]])) {
      transTypes[[p]] <- "upper-bounded"
      theta_t[,i] <- log(ub[[p]] - theta[,i])
    } else if (lb[[p]] < ub[[p]] && is.finite(lb[[p]]) && is.finite(ub[[p]])) {
      transTypes[[p]] <- "double-bounded"
      theta_t[,i] <- qnorm( (theta[,i] - lb[[p]])/(ub[[p]] - lb[[p]]) )
    } else {
      stop("Could not transform parameters, possibly due to invalid lower and/or upper
           prior bounds.")
    }

  }

  colnames(theta_t) <- paste0("trans_", colnames(theta))

  return(list(theta_t = theta_t, transTypes = transTypes))

}

#' @importFrom  stats pnorm
.invTransform2Real <- function(theta_t, lb, ub) {

  ### transform transformed samples back to original scales

  theta <- theta_t
  colnames(theta) <- stringr::str_sub(colnames(theta), 7)
  cn <- colnames(theta)

  for (i in seq_len(ncol(theta_t))) {

    p <- cn[i]

    if (lb[[p]] < ub[[p]] && is.infinite(lb[[p]]) && is.infinite(ub[[p]])) {
      theta[,i] <- theta_t[,i]
    } else if (lb[[p]] < ub[[p]] && is.finite(lb[[p]]) && is.infinite(ub[[p]])) {
      theta[,i] <- exp(theta_t[,i]) + lb[[p]]
    } else if (lb[[p]] < ub[[p]] && is.infinite(lb[[p]]) && is.finite(ub[[p]])) {
      theta[,i] <- ub[[p]] - exp(theta_t[,i])
    } else if (lb[[p]] < ub[[p]] && is.finite(lb[[p]]) && is.finite(ub[[p]])) {
      theta[,i] <- pnorm(theta_t[,i])*(ub[[p]] - lb[[p]]) + lb[[p]]
    } else {
      stop("Could not transform parameters, possibly due to invalid lower and/or upper
           prior bounds.")
    }

  }

  return(theta)

}

#' @importFrom  stats dnorm
.logJacobian <- function(theta_t, transTypes, lb, ub) {

  ### compute log of Jacobian

  logJ <- matrix(nrow = nrow(theta_t), ncol = ncol(theta_t))
  cn <- stringr::str_sub(colnames(theta_t), 7)

  for (i in seq_len( ncol(theta_t) )) {

    p <- cn[i]

    if (transTypes[[p]] == "unbounded") {
      logJ[,i] <- 0
    } else if (transTypes[[p]] == "lower-bounded") {
      logJ[,i] <- theta_t[,i]
    } else if (transTypes[[p]] == "upper-bounded") {
      logJ[,i] <- theta_t[,i]
    } else if (transTypes[[p]] == "double-bounded") {
      logJ[,i] <- log(ub[[p]] - lb[[p]]) + dnorm(theta_t[,i], log = TRUE)
    }

  }

  return(.rowSums(logJ, m = nrow(logJ), n = ncol(logJ)))

}

#' @importFrom stats median
.run.iterative.scheme <- function(q11, q12, q21, q22, r0, tol, L, method) {

  ### run iterative updating scheme (using "optimal" bridge function,
  ### see Meng & Wong, 1996)

  if (method == "normal") {
    l1 <- q11 - q12 # log(l)
    l2 <- q21 - q22 # log(ltilde)
  } else if (method == "warp3") {
    l1 <- -log(2) + log(det(L)) + (q11 - q12) # log(l)
    l2 <-  -log(2) + log(det(L)) + (q21 - q22) # log(ltilde)
  }

  lstar <- median(l1)
  n.1 <- length(l1)
  n.2 <- length(l2)
  s1 <- n.1/(n.1 + n.2)
  s2 <- n.2/(n.1 + n.2)
  rold <- -100
  r <- r0
  i <- 1

  e <- as.brob( exp(1) )

  while (abs((r - rold)/r) > tol) {

    cat(paste0("Iteration: ", i, "\n"))
    rold <- r
    numi <- as.numeric( e^(l2 - lstar)/(s1 * e^(l2 - lstar) + s2 *  r) )
    deni <- as.numeric( 1/(s1 * e^(l1 - lstar) + s2 * r) )
    r <- (n.1/n.2) * sum(numi)/sum(deni)
    i <- i + 1

  }

  logml <- log(r) + lstar
  return(list(logml = logml, niter = i-1))

}

#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom Matrix nearPD
#' @import Brobdingnag
#' @import stringr
#' @import coda
#' @import parallel
#' @import rlecuyer
#' @importFrom Rcpp sourceCpp
#' @import stats
.bridge.sampler.normal <- function(samples, log_posterior, data, lb, ub,
                                   cores, packages, varlist, rcppFile, r0, tol) {

  # transform parameters to real line
  tmp <- .transform2Real(samples, lb, ub)
  theta_t <- tmp$theta_t
  transTypes <- tmp$transTypes

  # split samples for proposal/iterative scheme
  nr <- nrow(samples)
  samples4fit_index <- seq_len(nr) %% 2 == TRUE # split samples in even/odd
  samples_4_fit <- theta_t[samples4fit_index, ,drop = FALSE]
  samples_4_iter <- theta_t[!samples4fit_index, , drop = FALSE]
  n_post <- nrow(samples_4_iter)

  # get mean & covariance matrix and generate samples from proposal
  m <- apply(samples_4_fit, 2, mean)
  V_tmp <- cov(samples_4_fit)
  V <- as.matrix(nearPD(V_tmp)$mat) # make sure that V is positive-definite
  gen_samples <- rmvnorm(n_post, mean = m, sigma = V)
  colnames(gen_samples) <- colnames(samples_4_iter)

  # evaluate multivariate normal distribution for posterior samples and generated samples
  q12 <- dmvnorm(samples_4_iter, mean = m, sigma = V, log = TRUE)
  q22 <- dmvnorm(gen_samples, mean = m, sigma = V, log = TRUE)

  # evaluate log likelihood times prior for posterior samples and generated samples
  if (cores == 1) {

    q11 <- apply(.invTransform2Real(samples_4_iter, lb, ub), 1, log_posterior,
                 data = data) + .logJacobian(samples_4_iter, transTypes, lb, ub)
    q21 <- apply(.invTransform2Real(gen_samples, lb, ub), 1, log_posterior,
                 data = data) + .logJacobian(gen_samples, transTypes, lb, ub)

  } else if (cores > 1) {

    cl <- parallel::makeCluster(cores, useXDR = FALSE)
    sapply(packages, function(x) parallel::clusterCall(cl = cl, "require", package = x, character.only = TRUE))
    parallel::clusterExport(cl = cl, varlist = varlist)

    if ( ! is.null(rcppFile)) {
      parallel::clusterExport(cl = cl, varlist = "rcppFile", envir = parent.frame())
      parallel::clusterCall(cl = cl, "require", package = "Rcpp", character.only = TRUE)
      parallel::clusterEvalQ(cl = cl, sourceCpp(file = rcppFile))
    }

    q11 <- parallel::parRapply(cl = cl, x = .invTransform2Real(samples_4_iter, lb, ub), log_posterior,
                               data = data) + .logJacobian(samples_4_iter, transTypes, lb, ub)
    q21 <- parallel::parRapply(cl = cl, x = .invTransform2Real(gen_samples, lb, ub), log_posterior,
                               data = data) + .logJacobian(gen_samples, transTypes, lb, ub)
    parallel::stopCluster(cl)

  }

  # run iterative updating scheme to compute log of marginal likelihood
  tmp <- .run.iterative.scheme(q11 = q11, q12 = q12, q21 = q21, q22 = q22,
                               r0 = r0, tol = tol, L = NULL, method = "normal")
  logml <- tmp$logml
  niter <- tmp$niter

  out <- list(logml = logml, niter = niter, method = "normal", q11 = q11,
              q12 = q12, q21 = q21, q22 = q22)

  return(out)

}

#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom Matrix nearPD
#' @import Brobdingnag
#' @import stringr
#' @import coda
#' @import rlecuyer
#' @importFrom Rcpp sourceCpp
#' @import stats
.bridge.sampler.warp3 <- function(samples, log_posterior, data, lb, ub,
                                  cores, packages, varlist, rcppFile, r0, tol) {

  # transform parameters to real line
  tmp <- .transform2Real(samples, lb, ub)
  theta_t <- tmp$theta_t
  transTypes <- tmp$transTypes

  # split samples for proposal/iterative scheme
  nr <- nrow(samples)
  samples4fit_index <- seq_len(nr) %% 2 == TRUE # split samples in even/odd
  samples_4_fit <- theta_t[samples4fit_index, ,drop = FALSE]
  samples_4_iter <- theta_t[!samples4fit_index, , drop = FALSE]
  n_post <- nrow(samples_4_iter)

  # get mean & covariance matrix and generate samples from proposal
  m <- apply(samples_4_fit, 2, mean)
  V_tmp <- cov(samples_4_fit)
  V <- as.matrix(nearPD(V_tmp)$mat) # make sure that V is positive-definite
  L <- t(chol(V))
  gen_samples <- rmvnorm(n_post, sigma = diag(ncol(samples_4_fit)))
  colnames(gen_samples) <- colnames(samples_4_iter)

  # evaluate multivariate normal distribution for posterior samples and generated samples
  q12 <- dmvnorm((samples_4_iter - matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE)) %*%
                   t(solve(L)), sigma = diag(ncol(samples_4_fit)), log = TRUE)
  q22 <- dmvnorm(gen_samples, sigma = diag(ncol(samples_4_fit)), log = TRUE)

  e <- as.brob( exp(1) )

  # evaluate log likelihood times prior for posterior samples and generated samples
  if (cores == 1) {

    q11 <- log(e^(apply(.invTransform2Real(samples_4_iter, lb, ub), 1, log_posterior,
                        data = data) + .logJacobian(samples_4_iter, transTypes, lb, ub)) +
                 e^(apply(.invTransform2Real(matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                               samples_4_iter, lb, ub), 1, log_posterior, data = data) +
                      .logJacobian(matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                     samples_4_iter, transTypes, lb, ub)))
    q21 <- log(e^(apply(.invTransform2Real(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                             gen_samples %*% t(L), lb, ub), 1, log_posterior, data = data) +
                    .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                   gen_samples %*% t(L), transTypes, lb, ub)) +
                 e^(apply(.invTransform2Real(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) +
                                               gen_samples %*% t(L), lb, ub), 1, log_posterior, data = data) +
                      .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) +
                                     gen_samples %*% t(L), transTypes, lb, ub)))

  } else if (cores > 1) {

    cl <- parallel::makeCluster(cores, useXDR = FALSE)
    sapply(packages, function(x) parallel::clusterCall(cl = cl, "require", package = x, character.only = TRUE))
    parallel::clusterExport(cl = cl, varlist = varlist)

    if ( ! is.null(rcppFile)) {
      parallel::clusterExport(cl = cl, varlist = "rcppFile", envir = parent.frame())
      parallel::clusterCall(cl = cl, "require", package = "Rcpp", character.only = TRUE)
      parallel::clusterEvalQ(cl = cl, sourceCpp(file = rcppFile))
    }

    q11 <- log(e^(parallel::parRapply(cl = cl, x = .invTransform2Real(samples_4_iter, lb, ub),
                                      log_posterior, data = data) +
                    .logJacobian(samples_4_iter, transTypes, lb, ub)) +
                 e^(parallel::parRapply(cl = cl, x = .invTransform2Real(matrix(2*m, nrow = n_post,
                                                                               ncol = length(m),
                                                                               byrow = TRUE) -
                                                                          samples_4_iter, lb, ub),
                                        log_posterior, data = data) +
                      .logJacobian(matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                     samples_4_iter, transTypes, lb, ub)))
    q21 <- log(e^(parallel::parRapply(cl = cl, x = .invTransform2Real(matrix(m, nrow = n_post,
                                                                             ncol = length(m), byrow = TRUE) -
                                               gen_samples %*% t(L), lb, ub), log_posterior, data = data) +
                    .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                   gen_samples %*% t(L), transTypes, lb, ub)) +
                 e^(parallel::parRapply(cl = cl, x = .invTransform2Real(matrix(m, nrow = n_post, ncol = length(m),
                                                                               byrow = TRUE) +
                                                 gen_samples %*% t(L), lb, ub), log_posterior, data = data) +
                      .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) +
                                     gen_samples %*% t(L), transTypes, lb, ub)))

    parallel::stopCluster(cl)

  }

  # run iterative updating scheme to compute log of marginal likelihood
  tmp <- .run.iterative.scheme(q11 = q11, q12 = q12, q21 = q21, q22 = q22,
                               r0 = r0, tol = tol, L = L, method = "warp3")
  logml <- tmp$logml
  niter <- tmp$niter

  out <- list(logml = logml, niter = niter, method = "warp3", q11 = q11,
              q12 = q12, q21 = q21, q22 = q22)

  return(out)

}

#' Computes log marginal likelihood via bridge sampling
#' @export
#' @title Computes log marginal likelihood via bridge sampling
#' @name bridge_sampler
#' @param samples matrix with posterior samples (colnames need to correspond to parameter names in \code{lb} and \code{ub}).
#' @param log_posterior function that takes a single row of \code{samples} and the \code{data} and returns the log of the unnormalized posterior density (i.e., a scalar value).
#' @param data data.
#' @param lb named vector with lower bounds for parameters.
#' @param ub named vector with upper bounds for parameters.
#' @param method either \code{"normal"} or \code{"warp3"}.
#' @param cores number of cores used for evaluating \code{log_posterior}.
#' @param packages character vector with names of packages needed for evaluating \code{log_posterior} in parallel (only relevant if \code{cores} > 1).
#' @param varlist character vector with names of variables needed for evaluating \code{log_posterior} (only needed if \code{cores} > 1 as those objects will be exported to the nodes). Those objects need to exist in the \code{\link{.GlobalEnv}}.
#' @param rcppFile in case \code{cores} > 1 and log_posterior is an Rcpp function, rcppFile specifies the path to the cpp file (needs to be compiled on all cores).
#' @details Bridge sampling is implemented as described in Meng and Wong (1996, see equation 4.1) using the "optimal" bridge function. When \code{method = "normal"}, the proposal distribution is a multivariate normal distribution with mean vector equal to the column means of \code{samples} and covariance matrix equal to the sample covariance matrix of \code{samples}. When \code{method = "warp3"},
#' the proposal distribution is a standard multivariate normal distribution and the posterior distribution is "warped" (Meng & Schilling, 2002) so that it has the same mean vector, covariance matrix and skew as samples. \code{method = "warp3"} usually yields more precise results, but it takes approximately twice as long as \code{method = "normal"}.
#' @return a list of class \code{"bridge"} with components:
#' \itemize{
#'  \item \code{logml}: estimate of log marginal likelihood.
#'  \item \code{niter}: number of iterations of the iterative updating scheme.
#'  \item \code{method}: bridge sampling method that was used to obtain estimate.
#'  \item \code{q11}: log_posterior evaluations for posterior samples.
#'  \item \code{q12}: log proposal evaluations for posterior samples.
#'  \item \code{q21}: log_posterior evaluations for samples from proposal.
#'  \item \code{q22}: log proposal evaluations for samples from proposal.
#' }
#' @author Quentin F. Gronau
#' @references
#' Meng, X.-L., & Wong, W. H. (1996). Simulating ratios of normalizing constants via a simple identity: A theoretical exploration. Statistica Sinica, 6, 831-860.
#'
#' Meng, X.-L., & Schilling, S. (2002). Warp bridge sampling. Journal of Computational and Graphical Statistics, 11(3), 552-586.
#'
#' @example examples/example.bridge_sampler.R
bridge_sampler <- function(samples = NULL, log_posterior = NULL, data = NULL,
                           lb = NULL, ub = NULL, method = "normal", cores = 1,
                           packages = NULL, varlist = NULL, rcppFile = NULL) {

  # see Meng & Wong (1996), equation 4.1

  out <- do.call(what = paste0(".bridge.sampler.", method),
                 args = list(samples = samples, log_posterior = log_posterior,
                             data = data, lb = lb, ub = ub, cores = cores,
                             packages = packages, varlist = varlist, rcppFile = rcppFile, r0 = 0,
                             tol = 1e-10))
  class(out) <- "bridge"
  return(out)

}

#' @export
print.bridge <- function(x, ...) {

  cat("Bridge sampling estimate of the log marginal likelihood: ",
      round(x$logml, 3), ". \nEstimate obtained in ", x$niter,
      " iterations via method \"", x$method, "\".", sep = "")
}


