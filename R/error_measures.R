#' Computes error measures for estimated marginal likelihood.
#' @export
#' @title Error Measures for Estimated Marginal Likelihood
#' @param bridge_object an object of class \code{"bridge"} or \code{"bridge_list"} as returned from \code{\link{bridge_sampler}}.
#' @param na.rm a logical indicating whether missing values in logml estimates should be removed.  Ignored for the \code{bridge} method.
#' @param ... additional arguments (currently ignored).
#' @details Computes error measures for marginal likelihood bridge sampling estimates. The approximate errors for a \code{bridge_object} of class \code{"bridge"} that has been obtained with \code{method = "normal"} and \code{repetitions = 1} are based on Fruehwirth-Schnatter (2004).
#' Not applicable in case the object of class \code{"bridge"} has been obtained with \code{method = "warp3"} and \code{repetitions = 1}.
#' To assess the uncertainty of the estimate in this case, it is recommended to run the \code{"warp3"} procedure multiple times.
#' @return If \code{bridge_object} is of class \code{"bridge"} and has been obtained with \code{method = "normal"} and \code{repetitions = 1}, returns a list with components:
#' \itemize{
#'  \item \code{re2}: approximate relative mean-squared error for marginal likelihood estimate.
#'  \item \code{cv}: approximate coefficient of variation for marginal likelihood estimate (assumes that bridge estimate is unbiased).
#'  \item \code{percentage}: approximate percentage error of marginal likelihood estimate.
#' }
#' If \code{bridge_object} is of class \code{"bridge_list"}, returns a list with components:
  #' \itemize{
  #'  \item \code{min}: minimum of the log marginal likelihood estimates.
  #'  \item \code{max}: maximum of the log marginal likelihood estimates.
  #'  \item \code{IQR}: interquartile range of the log marginal likelihood estimates.
  #' }
#' @author Quentin F. Gronau
#' @note For examples, see \code{\link{bridge_sampler}} and the accompanying vignettes: \cr \code{vignette("bridgesampling_example_jags")} \cr \code{vignette("bridgesampling_example_stan")}
#'
#' @seealso The \code{summary} methods for \code{bridge} and \code{bridge_list} objects automatically invoke this function, see \code{\link{bridge-methods}}.
#'
#' @references
#' Fruehwirth-Schnatter, S. (2004). Estimating marginal likelihoods for mixture and Markov switching models using bridge sampling techniques. \emph{The Econometrics Journal, 7}, 143-167. \url{http://dx.doi.org/10.1111/j.1368-423X.2004.00125.x}
#' @import Brobdingnag
#' @importFrom coda spectrum0.ar
#' @export
error_measures <- function (bridge_object, ...) {
  UseMethod("error_measures", bridge_object)
}

#' @rdname error_measures
#' @export
error_measures.bridge <- function(bridge_object,...) {

  if (bridge_object$method == "warp3") {
    stop(paste0("error_measures not implemented for warp3 method with",
                "\n  repetitions = 1.",
                "\n  We recommend to run the warp3 procedure multiple times",
                "\n  to assess the uncertainty of the estimate."))
  }

  e <- as.brob( exp(1) )

  ml <- e^(bridge_object$logml)
  g_p <- e^(bridge_object$q12)
  g_g <- e^(bridge_object$q22)
  priorTimesLik_p <- e^(bridge_object$q11)
  priorTimesLik_g <- e^(bridge_object$q21)
  p_p <- priorTimesLik_p/ml
  p_g <- priorTimesLik_g/ml

  N1 <- length(p_p)
  N2 <- length(g_g)
  s1 <- N1/(N1 + N2)
  s2 <- N2/(N1 + N2)

  f1 <- as.numeric( p_g/(s1*p_g + s2*g_g) )
  f2 <- as.numeric( g_p/(s1*p_p + s2*g_p) )
  rho_f2 <- spectrum0.ar( f2 )$spec

  term1 <- 1/N2 * var( f1 ) / mean( f1 )^2
  term2 <- rho_f2/N1 * var( f2 ) / mean( f2 )^2

  re2 <- term1 + term2

  # convert to coefficient of variation (assumes that bridge estimate is unbiased)
  cv <- sqrt(re2)

  # convert to percentage error
  percentage <- scales::percent(cv)
  return(list(re2 = re2, cv = cv, percentage = percentage))

}

#' @rdname error_measures
#' @export
error_measures.bridge_list <- function(bridge_object, na.rm = TRUE, ...) {

  return(list(min = min(bridge_object$logml, na.rm = na.rm),
              max = max(bridge_object$logml, na.rm = na.rm),
              IQR = stats::IQR(bridge_object$logml, na.rm = na.rm)))

}
