#' Computes posterior model probabilities.
#' @export
#' @title Posterior Model Probabilities from Marginal Likelihoods
#' @param ... objects of class \code{"bridge"} as returned from \code{\link{bridge_sampler}}.
#' @param prior_prob numeric vector with prior model probabilities. If omitted, a uniform prior is used (i.e., all models are equally likely a priori).
#' @return Numeric vector with posterior model probabilities.
#' @author Quentin F. Gronau
#' @note For examples, see \code{\link{bridge_sampler}} and the accompanying vignette: \cr \code{vignette("bridgesampling_example")}
compute_post_prob <- function(..., prior_prob) {

  dots <- list(...)
  mc <- match.call()
  if(!is.null(names(mc))) {
    mc <- mc[names(mc) != "prior_prob"]
  }

  mc <- mc[-1]

  if( ! all(vapply(dots, FUN = function(x) class(x) == "bridge", FUN.VALUE = TRUE)))
    stop("Not all passed objects are of class 'bridge'.")

  if(missing(prior_prob))
    prior_prob <- rep(1/length(dots), length(dots))

  if(!isTRUE(all.equal(sum(prior_prob), 1)))
    stop("Prior model probabilities do not sum to one.")

  if(length(mc) != length(prior_prob))
    stop("Number of objects given needs to match number of elements in prior_prob.")

  logml <- vapply(dots, function(x) x$logml, FUN.VALUE = 0)
  post_prob <- exp(logml)*prior_prob / sum(exp(logml)*prior_prob)
  names(post_prob) <- make.unique(as.character(mc))

  if(!isTRUE(all.equal(sum(post_prob), 1)))
    stop("Posterior model probabilities do not sum to one.")

  return(post_prob)

}
