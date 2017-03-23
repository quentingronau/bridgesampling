#' Computes error measures for estimated marginal likelihood.
#' @export
#' @title Computing error measures for estimated marginal likelihood
#' @param bridgeObject output from \code{\link{bridge_sampler}}.
#' @details Computes approximate error measures for marginal likelihood bridge sampling estimates. Based on Fruehwirth-Schnatter (2004).
#' @return a list with the objects:
#' \itemize{
#'  \item \code{re2}: approximate relative mean squared error for marginal likelihood estimate.
#'  \item \code{cv}: coefficient of variation for marginal likelihood estimate (assumes that bridge estimate is unbiased).
#'  \item \code{percentage}: percentage error of marginal likelihood estimate.
#' }
#' @author Quentin F. Gronau
#' @note For examples, see \code{\link{bridge_sampler}} and the accompanying vignette.
#' @references
#' Frühwirth‐Schnatter, S. (2004). Estimating marginal likelihoods for mixture and Markov switching models using bridge sampling techniques. The Econometrics Journal, 7, 143-167.
#' @import Brobdingnag
#' @importFrom coda spectrum0.ar
error_measures <- function(bridgeObject) {

  if (bridgeObject$method == "warp3")
    stop("error_measures not implemented for warp3 method. We recommend to run
         the warp3 procedure multiple times to assess the uncertainty of the
         estimate.")

  e <- as.brob( exp(1) )

  ml <- e^(bridgeObject$logml)
  g_p <- e^(bridgeObject$q12)
  g_g <- e^(bridgeObject$q22)
  priorTimesLik_p <- e^(bridgeObject$q11)
  priorTimesLik_g <- e^(bridgeObject$q21)
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
