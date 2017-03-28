#' Generic function that computes posterior model probabilities from marginal likelihoods
#' @export
#' @title Posterior Model Probabilities from Marginal Likelihoods
#' @param x Object of class \code{"bridge"} or \code{"bridge_list"} as returned from \code{\link{bridge_sampler}}. The default method simply assumes all passed objects are numeric log marginal likelihoods (e.g., from \code{\link{logml}}).
#' @param ... further objects of class \code{"bridge"} or \code{"bridge_list"} as returned from \code{\link{bridge_sampler}}. Or numeric values for the default method.
#' @param prior_prob numeric vector with prior model probabilities. If omitted, a uniform prior is used (i.e., all models are equally likely a priori).
#' @return Numeric vector with posterior model probabilities (with names derived from the input).
#' @author Quentin F. Gronau and Henrik Singmann
#' @note For realistic examples, see \code{\link{bridge_sampler}} and the accompanying vignette: \cr \code{vignette("bridgesampling_example")}
#' @example examples/example.post_prob.R
#'
post_prob <- function (x, ..., prior_prob) {
   UseMethod("post_prob", x)
}

#' @rdname post_prob
#' @export
post_prob.bridge <- function(x, ..., prior_prob) {
  dots <- list(...)
  mc <- match.call()
  modb <- (as.logical(vapply(dots, is, NA, "bridge")) |
             as.logical(vapply(dots, is, NA, "bridge_list"))
           )
  model_names <- c(deparse(mc[["x"]]), vapply(which(modb), function(x) deparse(mc[[x+2]]), ""))
  if (sum(modb) != length(dots))
    warning("Objects not of class 'bridge' or 'bridge_repetitions' are ignored.", call. = FALSE)

  logml <- vapply(c(list(x), dots[modb]), logml, FUN.VALUE = 0)

  .post_prob_calc(logml=logml, model_names = model_names, prior_prob=prior_prob)

}

#' @rdname post_prob
#' @export
post_prob.default <- function(x, ..., prior_prob) {
  dots <- list(...)
  mc <- match.call()
  model_names <- c(rep(deparse(mc[["x"]]), length(x)),
                   rep(vapply(seq_along(dots), function(x) deparse(mc[[x+2]]), ""),
                       times = vapply(dots, length, 0)))
  logml <- c(x, unlist(dots))
  if (!is.numeric(logml)) {
    stop("logml values need to be numeric", call. = FALSE)
  }
  .post_prob_calc(logml=logml, model_names = model_names, prior_prob=prior_prob)

}

.post_prob_calc <- function(logml, model_names, prior_prob) {
  e <- as.brob(exp(1))

  if(missing(prior_prob))
    prior_prob <- rep(1/length(logml), length(logml))

  if(!isTRUE(all.equal(sum(prior_prob), 1)))
    stop("Prior model probabilities do not sum to one.", call. = FALSE)

  if(length(logml) != length(prior_prob))
    stop("Number of objects/logml-values needs to match number of elements in prior_prob.", call. = FALSE)

  post_prob <- as.numeric(e^logml*prior_prob / sum(e^logml*prior_prob))
  names(post_prob) <- make.unique(as.character(model_names))

  if(!isTRUE(all.equal(sum(post_prob), 1)))
    warning("Posterior model probabilities do not sum to one.", call. = FALSE)

  return(post_prob)

}

