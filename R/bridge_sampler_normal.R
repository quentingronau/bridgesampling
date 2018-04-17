
.bridge.sampler.normal <- function(
  samples_4_fit, # matrix with already transformed samples for fitting the
                 # proposal (rows are samples), colnames are "trans_x" where
                 # x is the parameter name
  samples_4_iter, # matrix with already transformed samples for the
                  # iterative scheme (rows are samples), colnames are "trans_x"
                  # where x is the parameter name
  neff, # effective sample size of samples_4_iter (i.e., already transformed samples), scalar
  log_posterior,
  ...,
  data,
  lb, ub,
  transTypes, # types of transformations (unbounded/lower/upperbounded) for the different parameters (named character vector)
  param_types, # Sample space for transformations (real, circular, simplex)
  cores,
  repetitions,
  packages,
  varlist,
  envir,
  rcppFile,
  maxiter,
  silent,
  verbose,
  r0,
  tol1,
  tol2) {

  if (is.null(neff))
    neff <- nrow(samples_4_iter)

  n_post <- nrow(samples_4_iter)

  # get mean & covariance matrix and generate samples from proposal
  m <- apply(samples_4_fit, 2, mean)
  V_tmp <- cov(samples_4_fit)
  V <- as.matrix(nearPD(V_tmp)$mat) # make sure that V is positive-definite

  # sample from multivariate normal distribution and evaluate for posterior samples and generated samples
  q12 <- dmvnorm(samples_4_iter, mean = m, sigma = V, log = TRUE)
  gen_samples <- vector(mode = "list", length = repetitions)
  q22 <- vector(mode = "list", length = repetitions)
  for (i in seq_len(repetitions)) {
    gen_samples[[i]] <- rmvnorm(n_post, mean = m, sigma = V)
    colnames(gen_samples[[i]]) <- colnames(samples_4_iter)
    q22[[i]] <- dmvnorm(gen_samples[[i]], mean = m, sigma = V, log = TRUE)
  }

  # evaluate log of likelihood times prior for posterior samples and generated samples
  q21 <- vector(mode = "list", length = repetitions)
  if (cores == 1) {
    q11 <- apply(.invTransform2Real(samples_4_iter, lb, ub, param_types), 1, log_posterior,
                 data = data, ...) + .logJacobian(samples_4_iter, transTypes, lb, ub)
    for (i in seq_len(repetitions)) {
      q21[[i]] <- apply(.invTransform2Real(gen_samples[[i]], lb, ub, param_types), 1, log_posterior,
                        data = data, ...) + .logJacobian(gen_samples[[i]], transTypes, lb, ub)
    }
  } else if (cores > 1) {
    if ( .Platform$OS.type == "unix") {
      split1 <- .split_matrix(matrix=.invTransform2Real(samples_4_iter, lb, ub, param_types), cores=cores)
      q11 <- parallel::mclapply(split1, FUN =
                                  function(x) apply(x, 1, log_posterior, data = data, ...),
                                  mc.preschedule = FALSE,
                                  mc.cores = cores)
      q11 <- unlist(q11) + .logJacobian(samples_4_iter, transTypes, lb, ub)
      for (i in seq_len(repetitions)) {
        split2 <- .split_matrix(matrix=.invTransform2Real(gen_samples[[i]], lb, ub, param_types), cores = cores)
        q21[[i]] <- parallel::mclapply(split2, FUN =
                                  function(x) apply(x, 1, log_posterior, data = data, ...),
                                  mc.preschedule = FALSE,
                                  mc.cores = cores)
        q21[[i]] <- unlist(q21[[i]]) + .logJacobian(gen_samples[[i]], transTypes, lb, ub)
      }
    } else {
    cl <- parallel::makeCluster(cores, useXDR = FALSE)
    sapply(packages, function(x) parallel::clusterCall(cl = cl, "require", package = x,
                                                       character.only = TRUE))
    parallel::clusterExport(cl = cl, varlist = varlist, envir = envir)

    if ( ! is.null(rcppFile)) {
      parallel::clusterExport(cl = cl, varlist = "rcppFile", envir = parent.frame())
      parallel::clusterCall(cl = cl, "require", package = "Rcpp", character.only = TRUE)
      parallel::clusterEvalQ(cl = cl, Rcpp::sourceCpp(file = rcppFile))
    } else if (is.character(log_posterior)) {
      parallel::clusterExport(cl = cl, varlist = log_posterior, envir = envir)
    }

    q11 <- parallel::parRapply(cl = cl, x = .invTransform2Real(samples_4_iter, lb, ub, param_types), log_posterior,
                               data = data, ...) + .logJacobian(samples_4_iter, transTypes, lb, ub)
    for (i in seq_len(repetitions)) {
      q21[[i]] <- parallel::parRapply(cl = cl, x = .invTransform2Real(gen_samples[[i]], lb, ub, param_types), log_posterior,
                                      data = data, ...) + .logJacobian(gen_samples[[i]], transTypes, lb, ub)
    }
    parallel::stopCluster(cl)
    }
  }
  if (any(is.infinite(q11))) {
    warning(sum(is.infinite(q11)), " of the ", length(q11)," log_prob() evaluations on the posterior draws produced -Inf/Inf.", call. = FALSE)
  }
  for (i in seq_len(repetitions)) {
    if (any(is.infinite(q21[[i]]))) {
      warning(sum(is.infinite(q21[[i]])), " of the ", length(q21[[i]])," log_prob() evaluations on the proposal draws produced -Inf/Inf.", call. = FALSE)
    }
  }
  if (any(is.na(q11))) {
    warning(sum(is.na(q11)), " evaluation(s) of log_prob() on the posterior draws produced NA and have been replaced by -Inf.", call. = FALSE)
    q11[is.na(q11)] <- -Inf
  }
  for (i in seq_len(repetitions)) {
    if (any(is.na(q21[[i]]))) {
      warning(sum(is.na(q21[[i]])), " evaluation(s) of log_prob() on the proposal draws produced NA nd have been replaced by -Inf.", call. = FALSE)
      q21[[i]][is.na(q21[[i]])] <- -Inf
    }
  }
  if(verbose) {
    print("summary(q12): (log_dens of proposal for posterior samples)")
    print(summary(q12))
    print("summary(q22): (log_dens of proposal for generated samples)")
    print(lapply(q22, summary))
    print("summary(q11): (log_dens of posterior for posterior samples)")
    print(summary(q11))
    print("summary(q21): (log_dens of posterior for generated samples)")
    print(lapply(q21, summary))
  }

  logml <- numeric(repetitions)
  niter <- numeric(repetitions)
  # run iterative updating scheme to compute log of marginal likelihood
  for (i in seq_len(repetitions)) {
    tmp <- .run.iterative.scheme(q11 = q11, q12 = q12, q21 = q21[[i]], q22 = q22[[i]],
                                 r0 = r0, tol = tol1, L = NULL, method = "normal",
                                 maxiter = maxiter, silent = silent,
                                 criterion = "r", neff = neff)
    if (is.na(tmp$logml)) {
      warning("logml could not be estimated within maxiter, rerunning with adjusted starting value. \nEstimate might be more variable than usual.", call. = FALSE)
      lr <- length(tmp$r_vals)
      # use geometric mean as starting value
      r0_2 <- sqrt(tmp$r_vals[lr - 1]*tmp$r_vals[lr])
      tmp <- .run.iterative.scheme(q11 = q11, q12 = q12, q21 = q21[[i]], q22 = q22[[i]],
                                   r0 = r0_2, tol = tol2, L = NULL, method = "normal",
                                   maxiter = maxiter, silent = silent,
                                   criterion = "logml", neff = neff)
      tmp$niter <- maxiter + tmp$niter
    }

    logml[i] <- tmp$logml
    niter[i] <- tmp$niter
    if (niter[i] == maxiter)
      warning("logml could not be estimated within maxiter, returning NA.", call. = FALSE)
  }

  if (repetitions == 1) {
    out <- list(logml = logml, niter = niter, method = "normal", q11 = q11,
                q12 = q12, q21 = q21[[1]], q22 = q22[[1]])
    class(out) <- "bridge"
  } else if (repetitions > 1) {
    out <- list(logml = logml, niter = niter, method = "normal", repetitions = repetitions)
    class(out) <- "bridge_list"
  }

  return(out)

}
