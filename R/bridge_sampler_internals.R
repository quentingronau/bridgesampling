
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

    if (any(is.infinite(numi)) || any(is.infinite(deni))) {
      warning("Infinite value in iterative scheme, returning NA.\n Try rerunning with more samples.", call. = FALSE)
      return(list(logml = NA, niter = i))

    }

    r <- (n.1/n.2) * sum(numi)/sum(deni)
    i <- i + 1

  }

  logml <- log(r) + lstar

  if (i >= maxiter) {
    warning("logml could not be estimated within maxiter, returning NA.", call. = FALSE)
    return(list(logml = NA, niter = i-1))
  }


  return(list(logml = logml, niter = i-1))

}
