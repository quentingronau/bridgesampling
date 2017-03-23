#' Computes the Bayes factor for a model compairson between two models.
#' @export
#' @title Computing the Bayes factor
#' @param bridgeObject1 output from \code{\link{bridge_sampler}}.
#' @param bridgeObject2 output from \code{\link{bridge_sampler}}.
#' @param log Boolean. If \code{TRUE}, the function returns the log of the Bayes factor. Default is \code{FALSE}.
#' @details Computes the Bayes factor (Kass & Raftery, 1995) in favor of the model associated with \code{bridgeObject1} over the model associated with \code{bridgeObject2}.
#' @return Value of the Bayes factor (scalar).
#' @author Quentin F. Gronau
#' @note For examples, see \code{\link{bridge_sampler}} and the accompanying vignette: \code{vignette("bridgesampling_example")}.
#' @references
#' Kass, R. E., & Raftery, A. E. (1995). Bayes factors. Journal of the american statistical association, 90(430), 773-795.
compute_bf <- function(bridgeObject1, bridgeObject2, log = FALSE) {

  bf <- bridgeObject1$logml - bridgeObject2$logml

  if (! log)
    bf <- exp(bf)

  return(bf)

}
