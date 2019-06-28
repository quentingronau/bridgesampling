# Helper function to represent circular variables (such as mean directions) as
# "gapless" numerical representations.
.gaplessCircular <- function(th) {

  # Mean direction
  md     <- atan2(sum(sin(th)), sum(cos(th)))

  # Shift th so that it is unlikely to have a gap.
  ((th - md + pi) %% (2*pi)) - pi + md
}



#### for matrix method ######

.transform2Real <- function(theta, lb, ub,
                            theta_types = rep("real", ncol(theta))) {

  ### transform samples to real line

  theta_t <- theta
  transTypes <- character(ncol(theta))
  cn <- colnames(theta)
  names(theta_types) <- names(transTypes) <- cn

  # Because the simplex transform must be done on all simplex parameters at
  # once, do it before the loop. This transformation follows the Stan reference
  # manual. For simplex variables, we expect one parameter less than the number
  # of weights due to the contstraint sum(simplex_theta) == 1.
  is_simplex_theta <- theta_types == "simplex"
  if (any(is_simplex_theta)) {

    # Select the simplex variables
    simplex_theta <- theta[, is_simplex_theta, drop = FALSE]

    # Simplex dimensionality
    simdim <- ncol(simplex_theta)
    cs     <- cbind(0L, t(apply(simplex_theta, 1L, cumsum))[, -simdim, drop = FALSE])

    # Get the break proportions.
    z_k    <- (simplex_theta / (1L - cs))
    y_k    <- log(z_k) - log(1L - z_k) + matrix(log(simdim:1L),
                                                nrow(theta), simdim, byrow = TRUE)

    theta_t[, is_simplex_theta] <- y_k
    transTypes[is_simplex_theta] <- "simplex"
  }


  for (i in seq_len(ncol(theta))) {

    p <- cn[i]

    if (theta_types[[p]] == "circular") {
      transTypes[[p]] <- "circular"
      theta_t[,i] <- .gaplessCircular(theta[,i])

    } else if (theta_types[[p]] == "real") {
      if (any(theta[,i] < lb[[p]])) {
        stop("Parameter values (samples) cannot be smaller than lb: ", p,
             call. = FALSE)
      }
      if (any(theta[,i] > ub[[p]])) {
        stop("Parameter values (samples) cannot be larger than ub: ", p,
             call. = FALSE)
      }
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

        # Finally, give an error except for simplex variables, which are already
        # transformed.
      } else if (theta_types[p] != "simplex")
        stop(paste("Could not transform parameters, possibly due to invalid",
                   "lower and/or upper prior bounds."))
      }
    }


  colnames(theta_t) <- paste0("trans_", colnames(theta))

  return(list(theta_t = theta_t, transTypes = transTypes))
}

.invTransform2Real <- function(theta_t, lb, ub,
                               theta_types = rep("real", ncol(theta))) {

  ### transform transformed samples back to original scales

  theta <- theta_t
  colnames(theta) <- stringr::str_sub(colnames(theta), 7)
  cn <- colnames(theta)
  names(theta_types) <- cn

  # Because the simplex transform must be done on all simplex parameters at
  # once, do it before the loop. This transformation follows the Stan reference
  # manual. For simplex variables, we expect one parameter less than the number
  # of weights due to the contstraint sum(simplex_theta) == 1.
  is_simplex_theta <- theta_types == "simplex"
  if (any(is_simplex_theta)) {

    # Select the simplex variables
    simplex_theta <- theta_t[, is_simplex_theta, drop = FALSE]

    # Simplex dimensionality
    simdim <- ncol(simplex_theta)

    logitz <- simplex_theta - matrix(log(simdim:1L),
                                      nrow(theta), simdim, byrow = TRUE)
    z_k    <- exp(logitz) / (1 + exp(logitz))
    x_k <- z_k

    if (simdim > 1) {
      for (k in 2:simdim) {
        x_k[, k] <- (1 - rowSums(x_k[, 1:(k - 1), drop = FALSE])) * z_k[, k]
      }
    }

    theta[, is_simplex_theta] <- x_k
  }

  # Note that the circular variables are not transformed back, because they are
  # simply a different numerical representation.
  for (i in seq_len(ncol(theta_t))) {

    p <- cn[i]

    if (theta_types[[p]] == "real") {
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

  }

  return(theta)

}

.logJacobian <- function(theta_t, transTypes, lb, ub) {

  ### compute log of Jacobian

  logJ <- matrix(nrow = nrow(theta_t), ncol = ncol(theta_t))
  cn <- stringr::str_sub(colnames(theta_t), 7)

  # Separate the computations for the simplex
  is_simplex_theta <- transTypes == "simplex"
  if (any(is_simplex_theta)) {

    # Select the simplex variables
    simplex_theta <- theta_t[, is_simplex_theta, drop = FALSE]

    # Simplex dimensionality, this is K - 1
    simdim <- ncol(simplex_theta)

    logitz <- simplex_theta - matrix(log(simdim:1L),
                                     nrow(theta_t), simdim, byrow = TRUE)
    z_k <- exp(logitz) / (1 + exp(logitz))
    x_k <- z_k

    # Sum_x_k is the length of the remaining stick at step k. At the start, the
    # whole stick is still left
    sum_x_k <- matrix(nrow = nrow(theta_t), ncol = simdim)
    sum_x_k[, 1] <- 1

    if (simdim > 1) {
      for (k in 2:simdim) {
        rsx <- rowSums(x_k[, 1:(k - 1), drop = FALSE])
        x_k[, k]     <- (1 - rsx) * z_k[, k]
        sum_x_k[, k] <- (1 - rsx)
      }
    }
    logJ[, is_simplex_theta] <- log(z_k) + log(1 - z_k) + log(sum_x_k)
  }


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
    } else if (transTypes[[p]] == "circular") {
      logJ[,i] <- 0
    }
  }

  return(.rowSums(logJ, m = nrow(logJ), n = ncol(logJ)))

}

.split_matrix <- function(matrix, cores) {
  out <- vector("list", cores)
  borders <- ceiling(seq(from = 0, to = nrow(matrix), length.out = cores + 1))
  for (i in seq_len(cores)) {
    out[[i]] <- matrix[(borders[i] + 1):borders[i + 1], , drop = FALSE]
  }
  out
}

.run.iterative.scheme <- function(q11, q12, q21, q22, r0, tol, L,
                                  method, maxiter, silent,
                                  criterion, neff) {

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
  s1 <- neff/(neff + n.2)
  s2 <- n.2/(neff + n.2)
  r <- r0
  r_vals <- r
  logml <- log(r) + lstar
  logml_vals <-  logml
  criterion_val <- 1 + tol

  e <- as.brob( exp(1) )
  i <- 1

  while (i <= maxiter && criterion_val > tol) {

    if (! silent)
      cat(paste0("Iteration: ", i, "\n"))

    rold <- r
    logmlold <- logml
    numi <-  e^(l2 - lstar)/(s1 * e^(l2 - lstar) + s2 *  r)
    deni <- 1/(s1 * e^(l1 - lstar) + s2 * r)

    if (any(is.infinite(as.numeric(numi))) ||
        any(is.infinite(as.numeric((deni))))) {
      warning("Infinite value in iterative scheme, returning NA.\n Try rerunning with more samples.", call. = FALSE)
      return(list(logml = NA, niter = i))

    }

    r <- (n.1/n.2) * sum(numi)/sum(deni)
    r_vals <- c(r_vals, r)
    logml <- log(r) + lstar
    logml_vals <- c(logml_vals, logml)
    criterion_val <- switch(criterion, "r" = abs((r - rold)/r),
                            "logml" = abs((logml - logmlold)/logml))
    i <- i + 1

  }

  if (i >= maxiter) {
    return(list(logml = NA, niter = i-1, r_vals = r_vals))
  }

  return(list(logml = logml, niter = i-1))

}
