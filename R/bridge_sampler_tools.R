
#--------------------------------------------------------------------------
# functions for Stan support via rstan
#--------------------------------------------------------------------------

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
  out <- tryCatch(
    rstan::log_prob(object = data$stanfit, upars = s.row)
    , error = function(e) -Inf
  )
  if (is.na(out)) out <- -Inf
  return(out)
}

.restructure_upars <- function(upars, use_neff) {
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

  list(
    samples_4_fit = samples_4_fit
    , samples_4_iter = samples_4_iter
    , neff = neff
    , lb = lb, ub = ub
    , param_types = rep("real", ncol(samples_4_fit))
    , transTypes = transTypes
  )
}

.validate_cores <- function(cores) {
  # cores > 1 only for unix:
  if (!(.Platform$OS.type == "unix") & (cores != 1)) {
    warning("cores > 1 only possible on Unix/MacOs. Uses 'cores = 1' instead.", call. = FALSE)
    return(1L)
  } else {
    return(cores)
  }
}


#--------------------------------------------------------------------------
# functions for CmdStan support via cmdstanr
#--------------------------------------------------------------------------

.cmdstan_log_posterior <- function(s.row, data) {
  out <- tryCatch(
    data$stanfit$log_prob(unconstrained_variables = s.row)
    , error = function(e) -Inf
  )
  if (is.na(out)) out <- -Inf
  return(out)
}
