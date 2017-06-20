
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
  samples4fit_index <- seq_len(nr) %in% seq_len(round(nr/2)) # split samples in two parts
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
  if (any(is.infinite(q11))) {
    warning(sum(is.infinite(q11)), " of the ", length(q11)," log_prob() evaluations on the warp-transformed posterior draws produced -Inf/Inf.", call. = FALSE)
  }
  for (i in seq_len(repetitions)) {
    if (any(is.infinite(q21[[i]]))) {
      warning(sum(is.infinite(q21[[i]])), " of the ", length(q21[[i]])," log_prob() evaluations on the warp-transformed proposal draws produced -Inf/Inf.", call. = FALSE)
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
