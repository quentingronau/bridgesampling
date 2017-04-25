#' Generic function that computes Bayes factor(s) from marginal likelihoods
#' @export
#' @title Bayes Factor(s) from Marginal Likelihoods
#' @param x1 Object of class \code{"bridge"} or \code{"bridge_list"} as returned from \code{\link{bridge_sampler}}. Additionally, the default method assumes that \code{x1} is a single numeric log marginal likelihood (e.g., from \code{\link{logml}}) and will throw an error otherwise.
#' @param x2 Object of class \code{"bridge"} or \code{"bridge_list"} as returned from \code{\link{bridge_sampler}}. Additionally, the default method assumes that \code{x2} is a single numeric log marginal likelihood (e.g., from \code{\link{logml}}) and will throw an error otherwise.
#' @param log Boolean. If \code{TRUE}, the function returns the log of the Bayes factor. Default is \code{FALSE}.
#' @details Computes the Bayes factor (Kass & Raftery, 1995) in favor of the model associated with \code{x1} over the model associated with \code{x2}.
#' @return For the default method and the method for \code{"bridge"} objects, the (scalar) value of the Bayes factor in favor of the model associated with \code{x1}.
#'
#' For the method for \code{"bridge_list"} objects, a numeric vector consisting of Bayes factors where each element gives the Bayes factor for one set of logmls in favor of the model associated with \code{x1}. The length of this vector is given by the \code{"bridge_list"} element with the most \code{repetitions}. Elements with fewer repetitions will be recycled (with warning).
#' @author Quentin F. Gronau
#' @note For examples, see \code{\link{bridge_sampler}} and the accompanying vignette: \cr \code{vignette("bridgesampling_example")}
#' @references
#' Kass, R. E., & Raftery, A. E. (1995). Bayes factors. \emph{Journal of the American Statistical Association}, 90(430), 773-795. \url{http://dx.doi.org/10.1080/01621459.1995.10476572}
#' @importFrom methods is
bf <- function(x1, x2, log = FALSE) {
  UseMethod("bf", x1)
}

.bf_calc <- function(logml1, logml2, log) {
  bf <- logml1 - logml2
  if (! log)
    bf <- exp(bf)
  return(bf)
}

#' @rdname bf
#' @export
bf.bridge <- function(x1, x2, log = FALSE) {
  .bf_calc(logml(x1), logml(x2), log = log)
}


#' @rdname bf
#' @export
bf.bridge_list <- function(x1, x2, log = FALSE) {
  logml1 <- x1$logml
  logml2 <- x2$logml
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
bf.default <- function(x1, x2, log = FALSE) {
  if (!is.numeric(c(x1, x2))) {
    stop("logml values need to be numeric", call. = FALSE)
  }
  if (length(x1) > 1 || length(x2) > 1) {
    stop("Both logmls need to be scalar values.", call. = FALSE)
  }
  .bf_calc(x1, x2, log = log)
}
