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
#' @importFrom methods is
bf <- function(bridge_object1, bridge_object2, log = FALSE) {
  UseMethod("bf", bridge_object1)
}

.bf_calc <- function(logml1, logml2, log) {
  bf <- logml1 - logml2
  if (! log)
    bf <- exp(bf)
  return(bf)
}

#' @rdname bf
#' @export
bf.bridge <- function(bridge_object1, bridge_object2, log = FALSE) {
  .bf_calc(logml(bridge_object1), logml(bridge_object2), log = log)
}


#' @rdname bf
#' @export
bf.bridge_list <- function(bridge_object1, bridge_object2, log = FALSE) {
  logml1 <- bridge_object1$logml
  logml2 <- bridge_object2$logml
  len1 <- length(logml1)
  len2 <- length(logml2)
  max_len <- max(c(len1, len2))
  if (!all(c(len1, len2) == max_len)) {
    warning("Not all objects provide ", max_len, " logmls. Some values are recycled.", call. = FALSE)
    logml1 <- rep(logml1, length.out = max_len)
    logml2 <- rep(logml2, length.out = max_len)
  }
  .bf_calc(logml1, logml2, log = log)
}

#' @rdname bf
#' @export
bf.default <- function(bridge_object1, bridge_object2, log = FALSE) {
  logml1 <- bridge_object1
  logml2 <- bridge_object2
  .bf_calc(logml1, logml2, log = log)
}
