#' Generic function that computes Bayes factor(s) from marginal likelihoods. \code{bayes_factor()} is simply an (S3 generic) alias for \code{bf()}.
#' @export
#' @title Bayes Factor(s) from Marginal Likelihoods
#' @param x1 Object of class \code{"bridge"} or \code{"bridge_list"} as returned from \code{\link{bridge_sampler}}. Additionally, the default method assumes that \code{x1} is a single numeric log marginal likelihood (e.g., from \code{\link{logml}}) and will throw an error otherwise.
#' @param x2 Object of class \code{"bridge"} or \code{"bridge_list"} as returned from \code{\link{bridge_sampler}}. Additionally, the default method assumes that \code{x2} is a single numeric log marginal likelihood (e.g., from \code{\link{logml}}) and will throw an error otherwise.
#' @param log Boolean. If \code{TRUE}, the function returns the log of the Bayes factor. Default is \code{FALSE}.
#' @param ... currently not used here, but can be used by other methods.
#' @details Computes the Bayes factor (Kass & Raftery, 1995) in favor of the model associated with \code{x1} over the model associated with \code{x2}.
#' @return For the default method returns a list of class \code{"bf_default"} with components:
#' \itemize{
#'  \item \code{bf}: (scalar) value of the Bayes factor in favor of the model associated with \code{x1} over the model associated with \code{x2}.
#'  \item \code{log}: Boolean which indicates whether \code{bf} corresponds to the log Bayes factor.
#' }
#'
#'
#' For the method for \code{"bridge"} objects returns a list of class \code{"bf_bridge"} with components:
#' \itemize{
#'  \item \code{bf}:  (scalar) value of the Bayes factor in favor of the model associated with \code{x1} over the model associated with \code{x2}.
#'  \item \code{log}: Boolean which indicates whether \code{bf} corresponds to the log Bayes factor.
#' }
#'
#'
#' For the method for \code{"bridge_list"} objects returns a list of class \code{"bf_bridge_list"} with components:
#' \itemize{
#'  \item \code{bf}:  a numeric vector consisting of Bayes factors where each element gives the Bayes factor for one set of logmls in favor of the model associated with \code{x1} over the model associated with \code{x2}. The length of this vector is given by the \code{"bridge_list"} element with the most \code{repetitions}. Elements with fewer repetitions will be recycled (with warning).
#'  \item \code{bf_median_based}: (scalar) value of the Bayes factor in favor of the model associated with \code{x1} over the model associated with \code{x2} that is based on the median values of the logml estimates.
#'  \item \code{log}: Boolean which indicates whether \code{bf} corresponds to the log Bayes factor.
#' }
#' @author Quentin F. Gronau
#' @note For examples, see \code{\link{bridge_sampler}} and the accompanying vignettes: \cr \code{vignette("bridgesampling_example_jags")} \cr \code{vignette("bridgesampling_example_stan")}
#' @references
#' Kass, R. E., & Raftery, A. E. (1995). Bayes factors. \emph{Journal of the American Statistical Association, 90(430)}, 773-795. \url{http://dx.doi.org/10.1080/01621459.1995.10476572}
#' @importFrom methods is
bf <- function(x1, x2, log = FALSE, ...) {
  UseMethod("bf", x1)
}

#' @rdname bf
#' @export
bayes_factor <- function(x1, x2, log = FALSE, ...) {
  UseMethod("bayes_factor", x1)
}


#' @rdname bf
#' @export
bayes_factor.default <- function(x1, x2, log = FALSE, ...) {
  bf(x1, x2, log = log, ...)
}


.bf_calc <- function(logml1, logml2, log) {
  bf <- logml1 - logml2
  if (! log)
    bf <- exp(bf)
  return(bf)
}

#' @rdname bf
#' @export
bf.bridge <- function(x1, x2, log = FALSE, ...) {
  #name1 <- as.character(match.call()[[2]])
  #name2 <- as.character(match.call()[[3]])
  if (!inherits(x2, c("bridge", "bridge_list")))
    stop("x2 needs to be of class 'bridge' or 'bridge_list'.", call. = FALSE)
  bf <- .bf_calc(logml(x1), logml(x2), log = log)
  out <- list(bf = bf, log = log)
  class(out) <- "bf_bridge"
  return(out)
}


#' @rdname bf
#' @export
bf.bridge_list <- function(x1, x2, log = FALSE, ...) {
  if (!inherits(x2, c("bridge", "bridge_list")))
    stop("x2 needs to be of class 'bridge' or 'bridge_list'.", call. = FALSE)
  logml1 <- x1$logml
  logml2 <- x2$logml
  median1 <- median(logml1, na.rm = TRUE)
  median2 <- median(logml2, na.rm = TRUE)
  len1 <- length(logml1)
  len2 <- length(logml2)
  max_len <- max(c(len1, len2))
  if (!all(c(len1, len2) == max_len)) {
    warning("Not all objects provide ", max_len, " logmls. Some values are recycled.", call. = FALSE)
    logml1 <- rep(logml1, length.out = max_len)
    logml2 <- rep(logml2, length.out = max_len)
  }
  bf <- .bf_calc(logml1, logml2, log = log)
  bf_median_based <- .bf_calc(median1, median2, log = log)
  out <- list(bf = bf, bf_median_based = bf_median_based, log = log)
  class(out) <- "bf_bridge_list"
  return(out)
}

#' @rdname bf
#' @export
bf.default <- function(x1, x2, log = FALSE, ...) {
  if (!is.numeric(c(x1, x2))) {
    stop("logml values need to be numeric", call. = FALSE)
  }
  if (length(x1) > 1 || length(x2) > 1) {
    stop("Both logmls need to be scalar values.", call. = FALSE)
  }
  bf <- .bf_calc(x1, x2, log = log)
  out <- list(bf = bf, log = log)
  class(out) <- "bf_default"
  return(out)
}

######## Methods for bf objects:

#' @method print bf_bridge
#' @export
print.bf_bridge <- function(x, ...) {
  if (x$log) {
    cat("The estimated log Bayes factor in favor of x1 over x2 is equal to: ", round(x$bf, 5), "\n", sep = "")
  } else if (! x$log) {
    cat("The estimated Bayes factor in favor of x1 over x2 is equal to: ", round(x$bf, 5), "\n", sep = "")
  }
}

#' @method print bf_bridge_list
#' @export
print.bf_bridge_list <- function(x, na.rm = TRUE,...) {
  if (x$log) {
    cat("The estimated log Bayes factor in favor of x1 over x2 (based on the medians of the log marginal likelihood estimates) is equal to: ",
        round(x$bf_median_based, 5), "\nRange of estimates: ", round(range(x$bf, na.rm=na.rm)[1], 5), " to ",
        round(range(x$bf, na.rm = na.rm)[2], 5),
        "\nInterquartile range: ", round(stats::IQR(x$bf, na.rm = na.rm), 5), "\n", sep = "")
    if (any(is.na(x$bf))) warning(sum(is.na(x$bf))," log Bayes factor estimate(s) are NAs.", call. = FALSE)
  } else if (! x$log) {
    cat("The estimated Bayes factor in favor of x1 over x2 (based on the medians of the log marginal likelihood estimates) is equal to: ",
        round(x$bf_median_based, 5), "\nRange of estimates: ", round(range(x$bf, na.rm=na.rm)[1], 5), " to ",
        round(range(x$bf, na.rm = na.rm)[2], 5),
        "\nInterquartile range: ", round(stats::IQR(x$bf, na.rm = na.rm), 5), "\n", sep = "")
    if (any(is.na(x$bf))) warning(sum(is.na(x$bf))," Bayes factor estimate(s) are NAs.", call. = FALSE)
  }
}

#' @method print bf_default
#' @export
print.bf_default <- function(x, ...) {
  if (x$log) {
    cat("The log Bayes factor in favor of x1 over x2 is equal to: ", round(x$bf, 5), "\n", sep = "")
  } else if (! x$log) {
    cat("The Bayes factor in favor of x1 over x2 is equal to: ", round(x$bf, 5), "\n", sep = "")
  }
}
