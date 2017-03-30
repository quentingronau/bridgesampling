
#### for matrix method ######

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

.split_matrix <- function(matrix, cores) {
  out <- vector("list", cores)
  borders <- ceiling(seq(from = 0, to = nrow(matrix), length.out = cores+1))
  for (i in seq_len(cores)) {
    out[[i]] <- matrix[(borders[i]+1):borders[i+1],,drop = FALSE]
  }
  out
}
.run.iterative.scheme <- function(q11, q12, q21, q22, r0, tol, L,
                                  method, maxiter, silent) {

  ### run iterative updating scheme (using "optimal" bridge function,
  ### see Meng & Wong, 1996)

  if (method == "normal") {
    l1 <- q11 - q12 # log(l)
    l2 <- q21 - q22 # log(ltilde)
  } else if (method == "warp3") {
    l1 <- -log(2) + determinant(L)$modulus + (q11 - q12) # log(l)
    l2 <-  -log(2) + determinant(L)$modulus + (q21 - q22) # log(ltilde)
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

  while (i <= maxiter && abs((r - rold)/r) > tol) {

    if (! silent)
      cat(paste0("Iteration: ", i, "\n"))

    rold <- r
    numi <- as.numeric( e^(l2 - lstar)/(s1 * e^(l2 - lstar) + s2 *  r) )
    deni <- as.numeric( 1/(s1 * e^(l1 - lstar) + s2 * r) )

    if (any(is.infinite(numi)) || any(is.infinite(deni)))
      stop("Infinite value in iterative scheme. Try rerunning with more samples.", call. = FALSE)

    r <- (n.1/n.2) * sum(numi)/sum(deni)
    i <- i + 1

  }

  logml <- log(r) + lstar

  if (i >= maxiter)
    stop("logml could not be estimated within maxiter")

  return(list(logml = logml, niter = i-1))

}

.bridge.sampler.normal <- function(samples, log_posterior, ..., data, lb, ub,
                                   cores, repetitions, packages, varlist, envir,
                                   rcppFile, maxiter, silent, verbose,
                                   r0, tol) {

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
    q11 <- apply(.invTransform2Real(samples_4_iter, lb, ub), 1, log_posterior,
                 data = data, ...) + .logJacobian(samples_4_iter, transTypes, lb, ub)
    for (i in seq_len(repetitions)) {
      q21[[i]] <- apply(.invTransform2Real(gen_samples[[i]], lb, ub), 1, log_posterior,
                        data = data, ...) + .logJacobian(gen_samples[[i]], transTypes, lb, ub)
    }
  } else if (cores > 1) {
    if ( .Platform$OS.type == "unix") {
      split1 <- .split_matrix(matrix=.invTransform2Real(samples_4_iter, lb, ub), cores=cores)
      q11 <- parallel::mclapply(split1, FUN =
                                  function(x) apply(x, 1, log_posterior, data = data, ...),
                                  mc.preschedule = FALSE,
                                  mc.cores = cores)
      q11 <- unlist(q11) + .logJacobian(samples_4_iter, transTypes, lb, ub)
      for (i in seq_len(repetitions)) {
        split2 <- .split_matrix(matrix=.invTransform2Real(gen_samples[[i]], lb, ub), cores = cores)
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

    q11 <- parallel::parRapply(cl = cl, x = .invTransform2Real(samples_4_iter, lb, ub), log_posterior,
                               data = data, ...) + .logJacobian(samples_4_iter, transTypes, lb, ub)
    for (i in seq_len(repetitions)) {
      q21[[i]] <- parallel::parRapply(cl = cl, x = .invTransform2Real(gen_samples[[i]], lb, ub), log_posterior,
                                      data = data, ...) + .logJacobian(gen_samples[[i]], transTypes, lb, ub)
    }
    parallel::stopCluster(cl)
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
                                 r0 = r0, tol = tol, L = NULL, method = "normal",
                                 maxiter = maxiter, silent = silent)
    logml[i] <- tmp$logml
    niter[i] <- tmp$niter
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

.bridge.sampler.warp3 <- function(samples, log_posterior, ..., data, lb, ub,
                                  repetitions, cores, packages, varlist, envir,
                                  rcppFile, maxiter, silent, verbose,
                                  r0, tol) {

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

  # sample from multivariate normal distribution and evaluate for posterior samples and generated samples
  q12 <- dmvnorm((samples_4_iter - matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE)) %*%
                   t(solve(L)), sigma = diag(ncol(samples_4_fit)), log = TRUE)
  q22 <- vector(mode = "list", length = repetitions)
  gen_samples <- vector(mode = "list", length = repetitions)
  for (i in seq_len(repetitions)) {
    gen_samples[[i]] <- rmvnorm(n_post, sigma = diag(ncol(samples_4_fit)))
    colnames(gen_samples[[i]]) <- colnames(samples_4_iter)
    q22[[i]] <- dmvnorm(gen_samples[[i]], sigma = diag(ncol(samples_4_fit)), log = TRUE)
  }

  e <- as.brob( exp(1) )

  # evaluate log of likelihood times prior for posterior samples and generated samples
  q21 <- vector(mode = "list", length = repetitions)
  if (cores == 1) {

    q11 <- log(e^(apply(.invTransform2Real(samples_4_iter, lb, ub), 1, log_posterior,
                        data = data,...) + .logJacobian(samples_4_iter, transTypes, lb, ub)) +
                 e^(apply(.invTransform2Real(matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                               samples_4_iter, lb, ub), 1, log_posterior, data = data, ...) +
                      .logJacobian(matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                     samples_4_iter, transTypes, lb, ub)))
    for (i in seq_len(repetitions)) {
      q21[[i]] <- log(e^(apply(.invTransform2Real(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                             gen_samples[[i]] %*% t(L), lb, ub), 1, log_posterior, data = data, ...) +
                    .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                   gen_samples[[i]] %*% t(L), transTypes, lb, ub)) +
                 e^(apply(.invTransform2Real(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) +
                                               gen_samples[[i]] %*% t(L), lb, ub), 1, log_posterior, data = data, ...) +
                      .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) +
                                     gen_samples[[i]] %*% t(L), transTypes, lb, ub)))
    }

  } else if (cores > 1) {
        if ( .Platform$OS.type == "unix") {
      split1a <- .split_matrix(matrix=.invTransform2Real(samples_4_iter, lb, ub), cores=cores)
      split1b <- .split_matrix(matrix=.invTransform2Real(
        matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) - samples_4_iter, lb, ub
        ), cores=cores)
      q11a <- parallel::mclapply(split1a, FUN =
                                   function(x) apply(x, 1, log_posterior, data = data, ...),
                                 mc.preschedule = FALSE,
                                 mc.cores = cores)
      q11b <- parallel::mclapply(split1b, FUN =
                                   function(x) apply(x, 1, log_posterior, data = data, ...),
                                 mc.preschedule = FALSE,
                                 mc.cores = cores)
      q11 <- log(e^(unlist(q11a) + .logJacobian(samples_4_iter, transTypes, lb, ub)) +
                   e^(unlist(q11b) +
                        .logJacobian(matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                       samples_4_iter, transTypes, lb, ub)))

      for (i in seq_len(repetitions)) {
        tmp_mat2 <-  gen_samples[[i]] %*% t(L)
        split2a <- .split_matrix(matrix=.invTransform2Real(
          matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) - tmp_mat2, lb, ub
        ), cores=cores)
        split2b <- .split_matrix(matrix=.invTransform2Real(
          matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) + tmp_mat2, lb, ub
        ), cores=cores)
        q21a <- parallel::mclapply(split2a, FUN =
                                     function(x) apply(x, 1, log_posterior, data = data, ...),
                                   mc.preschedule = FALSE,
                                   mc.cores = cores)
        q21b <- parallel::mclapply(split2b, FUN =
                                     function(x) apply(x, 1, log_posterior, data = data, ...),
                                   mc.preschedule = FALSE,
                                   mc.cores = cores)
        q21[[i]] <- log(e^(unlist(q21a) +
                        .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                       tmp_mat2, transTypes, lb, ub)) +
                     e^(unlist(q21b) +
                          .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) +
                                         tmp_mat2, transTypes, lb, ub)))
      }
        } else {
          cl <- parallel::makeCluster(cores, useXDR = FALSE)
          sapply(packages, function(x) parallel::clusterCall(cl = cl, "require", package = x, character.only = TRUE))

          parallel::clusterExport(cl = cl, varlist = varlist, envir = envir)

          if ( ! is.null(rcppFile)) {
            parallel::clusterExport(cl = cl, varlist = "rcppFile", envir = parent.frame())
            parallel::clusterCall(cl = cl, "require", package = "Rcpp", character.only = TRUE)
            parallel::clusterEvalQ(cl = cl, Rcpp::sourceCpp(file = rcppFile))
          } else if (is.character(log_posterior)) {
            parallel::clusterExport(cl = cl, varlist = log_posterior, envir = envir)
          }

          q11 <- log(e^(parallel::parRapply(cl = cl, x = .invTransform2Real(samples_4_iter, lb, ub),
                                            log_posterior, data = data, ...) +
                          .logJacobian(samples_4_iter, transTypes, lb, ub)) +
                       e^(parallel::parRapply(cl = cl,
                                              x = .invTransform2Real(matrix(2*m, nrow = n_post,
                                                                            ncol = length(m), byrow = TRUE) -
                                                                       samples_4_iter, lb, ub),
                                              log_posterior, data = data, ...) +
                            .logJacobian(matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                           samples_4_iter, transTypes, lb, ub)))
          for (i in seq_len(repetitions)) {
            q21[[i]] <- log(e^(parallel::parRapply(cl = cl,
                                            x = .invTransform2Real(matrix(m, nrow = n_post,
                                                                          ncol = length(m), byrow = TRUE) -
                                                                     gen_samples[[i]] %*% t(L), lb, ub),
                                            log_posterior, data = data, ...) +
                          .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) -
                                         gen_samples[[i]] %*% t(L), transTypes, lb, ub)) +
                       e^(parallel::parRapply(cl = cl,
                                              x = .invTransform2Real(matrix(m, nrow = n_post,
                                                                            ncol = length(m),byrow = TRUE) +
                                                                       gen_samples[[i]] %*% t(L), lb, ub),
                                              log_posterior, data = data, ...) +
                            .logJacobian(matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) +
                                           gen_samples[[i]] %*% t(L), transTypes, lb, ub)))
          }
          parallel::stopCluster(cl)
        }

  }
  if (any(is.na(q11))) {
    warning(sum(is.na(q11)), " evaluation(s) of log_prob() on the warp-transformed posterior draws produced NA and have been replaced by -Inf.", call. = FALSE)
    q11[is.na(q11)] <- -Inf
  }
  for (i in seq_len(repetitions)) {
    if (any(is.na(q21[[i]]))) {
      warning(sum(is.na(q21[[i]])), " evaluation(s) of log_prob() on the warp-transformed proposal draws produced NA nd have been replaced by -Inf.", call. = FALSE)
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
                                 r0 = r0, tol = tol, L = L, method = "warp3",
                                 maxiter = maxiter, silent = silent)
    logml[i] <- tmp$logml
    niter[i] <- tmp$niter
  }
  if (repetitions == 1) {
    out <- list(logml = logml, niter = niter, method = "warp3", q11 = q11,
                q12 = q12, q21 = q21[[1]], q22 = q22[[1]])
    class(out) <- "bridge"
  } else if (repetitions > 1) {
    out <- list(logml = logml, niter = niter, method = "warp3", repetitions = repetitions)
    class(out) <- "bridge_list"
  }

  return(out)

}

######### for stanfit method ########

# taken from rstan:
.rstan_relist <- function (x, skeleton) {
  lst <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) dim(lst[[i]]) <- dim(skeleton[[i]])
  lst
}

# taken from rstan:
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
  tryCatch(rstan::log_prob(object = data$stanfit, upars = s.row), error = function(e) NA)
}
