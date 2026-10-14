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
  lb,
  ub,
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
  use_ess,
  r0,
  tol1,
  tol2,
  gradients_4_fit = NULL, # optional n_fit x p matrix of posterior scores
  # grad log p(theta) at the rows of samples_4_fit
  proposal_fit = c("sample", "hybrid")
) {
  proposal_fit <- match.arg(proposal_fit)

  if (is.null(neff)) {
    neff <- nrow(samples_4_iter)
  }

  n_post <- nrow(samples_4_iter)

  # get mean & covariance matrix and generate samples from proposal
  m <- apply(samples_4_fit, 2, mean)
  V_sample <- cov(samples_4_fit)

  use_gradients <- !is.null(gradients_4_fit) && proposal_fit != "sample"
  if (use_gradients) {
    if (
      nrow(gradients_4_fit) != nrow(samples_4_fit) ||
        ncol(gradients_4_fit) != ncol(samples_4_fit)
    ) {
      warning(
        "gradients_4_fit dimensions do not match samples_4_fit; ",
        "falling back to the sample covariance.",
        call. = FALSE
      )
      use_gradients <- FALSE
    }
  }

  # combiner records how V_sample and V_score were combined:
  #   "none"      : V = V_sample
  #   "geometric" : V = V_sample # V_score
  #   "geometric_fallback_arithmetic" : .geometric_mean_psd() fell back
  proposal_fit_info <- list(
    proposal_fit = proposal_fit,
    combiner = NA_character_,
    n_gradients_used = NA_integer_
  )

  if (use_gradients) {
    # Score-matching estimate of the inverse covariance. The scores are
    # centred: the uncentred form retains an O(bar s bar s^T) bias in
    # finite samples even though E_pi[s] = 0.
    finite_rows <- apply(is.finite(gradients_4_fit), 1, all)
    if (sum(finite_rows) < ncol(gradients_4_fit) + 1L) {
      warning(
        sprintf(
          "only %d of %d gradient rows are finite (need more than ncol = %d); ",
          sum(finite_rows),
          nrow(gradients_4_fit),
          ncol(gradients_4_fit)
        ),
        "falling back to the sample covariance.",
        call. = FALSE
      )
      use_gradients <- FALSE
    } else {
      if (sum(finite_rows) < nrow(gradients_4_fit)) {
        warning(
          sprintf(
            "%d of %d draws produced non-finite gradients; ",
            sum(!finite_rows),
            nrow(gradients_4_fit)
          ),
          "they are dropped from the score-matched proposal fit.",
          call. = FALSE
        )
      }
      G <- gradients_4_fit[finite_rows, , drop = FALSE]
      G <- sweep(G, 2L, colMeans(G), check.margin = FALSE)
      Prec_score <- crossprod(G) / nrow(G)
      Prec_score <- as.matrix(nearPD(Prec_score)$mat)
      V_score <- tryCatch(
        chol2inv(chol(Prec_score)),
        error = function(e) solve(Prec_score)
      )
      proposal_fit_info$n_gradients_used <- nrow(G)
    }
  }

  if (use_gradients) {
    # .geometric_mean_psd() expects positive-definite inputs; V_score is
    # already nearPD'd above, V_sample is not.
    V_sample_pd <- as.matrix(nearPD(V_sample)$mat)
    V_tmp <- .geometric_mean_psd(V_sample_pd, V_score)
    proposal_fit_info$combiner <- if (isTRUE(attr(V_tmp, "fallback"))) {
      "geometric_fallback_arithmetic"
    } else {
      "geometric"
    }
    attr(V_tmp, "fallback") <- NULL
    if (verbose) {
      cat(sprintf(
        "[proposal_fit = hybrid] combiner = %s\n",
        proposal_fit_info$combiner
      ))
      cat(sprintf(
        "  sample Sigma eigenvalues: %.3g .. %.3g\n",
        min(eigen(V_sample, symmetric = TRUE, only.values = TRUE)$values),
        max(eigen(V_sample, symmetric = TRUE, only.values = TRUE)$values)
      ))
      cat(sprintf(
        "  score Sigma eigenvalues : %.3g .. %.3g\n",
        min(eigen(V_score, symmetric = TRUE, only.values = TRUE)$values),
        max(eigen(V_score, symmetric = TRUE, only.values = TRUE)$values)
      ))
    }
  } else {
    V_tmp <- V_sample
    proposal_fit_info$combiner <- "none"
  }

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
    q11 <- apply(
      .invTransform2Real(samples_4_iter, lb, ub, param_types),
      1,
      log_posterior,
      data = data,
      ...
    ) +
      .logJacobian(samples_4_iter, transTypes, lb, ub)
    for (i in seq_len(repetitions)) {
      q21[[i]] <- apply(
        .invTransform2Real(gen_samples[[i]], lb, ub, param_types),
        1,
        log_posterior,
        data = data,
        ...
      ) +
        .logJacobian(gen_samples[[i]], transTypes, lb, ub)
    }
  } else if (cores > 1) {
    if (.Platform$OS.type == "unix") {
      split1 <- .split_matrix(
        matrix = .invTransform2Real(samples_4_iter, lb, ub, param_types),
        cores = cores
      )
      q11 <- parallel::mclapply(
        split1,
        FUN = function(x) apply(x, 1, log_posterior, data = data, ...),
        mc.preschedule = FALSE,
        mc.cores = cores
      )
      q11 <- unlist(q11) + .logJacobian(samples_4_iter, transTypes, lb, ub)
      for (i in seq_len(repetitions)) {
        split2 <- .split_matrix(
          matrix = .invTransform2Real(gen_samples[[i]], lb, ub, param_types),
          cores = cores
        )
        q21[[i]] <- parallel::mclapply(
          split2,
          FUN = function(x) apply(x, 1, log_posterior, data = data, ...),
          mc.preschedule = FALSE,
          mc.cores = cores
        )
        q21[[i]] <- unlist(q21[[i]]) +
          .logJacobian(gen_samples[[i]], transTypes, lb, ub)
      }
    } else {
      cl <- parallel::makeCluster(cores, useXDR = FALSE)
      sapply(packages, function(x) {
        parallel::clusterCall(
          cl = cl,
          "require",
          package = x,
          character.only = TRUE
        )
      })
      parallel::clusterExport(cl = cl, varlist = varlist, envir = envir)

      if (!is.null(rcppFile)) {
        parallel::clusterExport(
          cl = cl,
          varlist = "rcppFile",
          envir = parent.frame()
        )
        parallel::clusterCall(
          cl = cl,
          "require",
          package = "Rcpp",
          character.only = TRUE
        )
        parallel::clusterEvalQ(cl = cl, Rcpp::sourceCpp(file = rcppFile))
      } else if (is.character(log_posterior)) {
        parallel::clusterExport(cl = cl, varlist = log_posterior, envir = envir)
      }

      q11 <- parallel::parRapply(
        cl = cl,
        x = .invTransform2Real(samples_4_iter, lb, ub, param_types),
        log_posterior,
        data = data,
        ...
      ) +
        .logJacobian(samples_4_iter, transTypes, lb, ub)
      for (i in seq_len(repetitions)) {
        q21[[i]] <- parallel::parRapply(
          cl = cl,
          x = .invTransform2Real(gen_samples[[i]], lb, ub, param_types),
          log_posterior,
          data = data,
          ...
        ) +
          .logJacobian(gen_samples[[i]], transTypes, lb, ub)
      }
      parallel::stopCluster(cl)
    }
  }
  if (verbose) {
    print(
      "summary(q12): (log_dens of proposal (i.e., with dmvnorm) for posterior samples)"
    )
    print(summary(q12))
    print(
      "summary(q22): (log_dens of proposal (i.e., with dmvnorm) for generated samples)"
    )
    print(lapply(q22, summary))
    print(
      "summary(q11): (log_dens of posterior (i.e., with log_posterior) for posterior samples)"
    )
    print(summary(q11))
    print(
      "summary(q21): (log_dens of posterior (i.e., with log_posterior) for generated samples)"
    )
    print(lapply(q21, summary))
    .PROPOSALS <- vector("list", repetitions)
    # for (i in seq_len(repetitions)) {
    #   .PROPOSALS[[i]] <- .invTransform2Real(gen_samples[[i]], lb, ub, param_types)
    # }
    # assign(".PROPOSALS", .PROPOSALS, pos = .GlobalEnv)
    # message("All proposal samples written to .GlobalEnv as .PROPOSALS")
  }
  if (any(is.infinite(q11))) {
    warning(
      sum(is.infinite(q11)),
      " of the ",
      length(q11),
      " log_prob() evaluations on the posterior draws produced -Inf/Inf.",
      call. = FALSE
    )
  }
  for (i in seq_len(repetitions)) {
    if (any(is.infinite(q21[[i]]))) {
      warning(
        sum(is.infinite(q21[[i]])),
        " of the ",
        length(q21[[i]]),
        " log_prob() evaluations on the proposal draws produced -Inf/Inf.",
        call. = FALSE
      )
    }
  }
  if (any(is.na(q11))) {
    warning(
      sum(is.na(q11)),
      " evaluation(s) of log_prob() on the posterior draws produced NA and have been replaced by -Inf.",
      call. = FALSE
    )
    q11[is.na(q11)] <- -Inf
  }
  for (i in seq_len(repetitions)) {
    if (all(is.na(q21[[i]]))) {
      stop(
        "Evaluations of log_prob() on all proposal draws produced NA.\n",
        "E.g., rounded to 3 digits (use verbose = TRUE for all proposal samples):\n",
        deparse(
          round(
            .invTransform2Real(gen_samples[[i]], lb, ub, param_types)[1, ],
            3
          ),
          width.cutoff = 500L
        ),
        call. = FALSE
      )
    }
    if (any(is.na(q21[[i]]))) {
      warning(
        sum(is.na(q21[[i]])),
        " evaluation(s) of log_prob() on the proposal draws produced NA and have been replaced by -Inf.",
        call. = FALSE
      )
      q21[[i]][is.na(q21[[i]])] <- -Inf
    }
  }
  logml <- numeric(repetitions)
  niter <- numeric(repetitions)
  mcse_logmls <- numeric(repetitions)
  # run iterative updating scheme to compute log of marginal likelihood
  for (i in seq_len(repetitions)) {
    tmp <- .run.iterative.scheme(
      q11 = q11,
      q12 = q12,
      q21 = q21[[i]],
      q22 = q22[[i]],
      r0 = r0,
      tol = tol1,
      L = NULL,
      method = "normal",
      maxiter = maxiter,
      silent = silent,
      use_ess = use_ess,
      criterion = "r",
      neff = neff
    )
    if (is.na(tmp$logml) & !is.null(tmp$r_vals)) {
      warning(
        "logml could not be estimated within maxiter, rerunning with adjusted starting value. \nEstimate might be more variable than usual.",
        call. = FALSE
      )
      lr <- length(tmp$r_vals)
      # use geometric mean as starting value
      r0_2 <- sqrt(tmp$r_vals[[lr - 1]] * tmp$r_vals[[lr]])
      tmp <- .run.iterative.scheme(
        q11 = q11,
        q12 = q12,
        q21 = q21[[i]],
        q22 = q22[[i]],
        r0 = r0_2,
        tol = tol2,
        L = NULL,
        method = "normal",
        maxiter = maxiter,
        silent = silent,
        use_ess = use_ess,
        criterion = "logml",
        neff = neff
      )
      tmp$niter <- maxiter + tmp$niter
    }

    logml[i] <- tmp$logml
    mcse_logmls[i] <- tmp$mcse_logml
    niter[i] <- tmp$niter
    if (niter[i] == maxiter) {
      warning(
        "logml could not be estimated within maxiter, returning NA.",
        call. = FALSE
      )
    }
  }

  if (repetitions == 1) {
    out <- list(
      logml = logml,
      niter = niter,
      method = "normal",
      q11 = q11,
      q12 = q12,
      q21 = q21[[1]],
      q22 = q22[[1]],
      mcse_logml = mcse_logmls,
      proposal_fit_info = proposal_fit_info
    )
    class(out) <- "bridge"
  } else if (repetitions > 1) {
    out <- list(
      logml = logml,
      niter = niter,
      method = "normal",
      repetitions = repetitions,
      mcse_logml = mcse_logmls,
      proposal_fit_info = proposal_fit_info
    )
    class(out) <- "bridge_list"
  }

  return(out)
}
