#' Computes the Bayes factor for a model comparison between two models.
#' @export
#' @title Bayes Factor from Marginal Likelihoods
#' @param bridge_object1 an object of class \code{"bridge"} as returned from \code{\link{bridge_sampler}}.
#' @param bridge_object2 an object of class \code{"bridge"} as returned from \code{\link{bridge_sampler}}.
#' @param log Boolean. If \code{TRUE}, the function returns the log of the Bayes factor. Default is \code{FALSE}.
#' @details Computes the Bayes factor (Kass & Raftery, 1995) in favor of the model associated with \code{bridge_object1} over the model associated with \code{bridge_object2}.
#' @return Value of the Bayes factor (scalar).
#' @author Quentin F. Gronau
#' @note For examples, see \code{\link{bridge_sampler}} and the accompanying vignette: \cr \code{vignette("bridgesampling_example")}
#' @references
#' Kass, R. E., & Raftery, A. E. (1995). Bayes factors. \emph{Journal of the American Statistical Association}, 90(430), 773-795. \url{http://dx.doi.org/10.1080/01621459.1995.10476572}
compute_bf <- function(bridge_object1, bridge_object2, log = FALSE) {

  bf <- bridge_object1$logml - bridge_object2$logml

  if (! log)
    bf <- exp(bf)

  return(bf)

}
