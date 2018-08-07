#' Generic function that computes posterior model probabilities from marginal
#' likelihoods.
#' @export
#' @title Posterior Model Probabilities from Marginal Likelihoods
#' @param x Object of class \code{"bridge"} or \code{"bridge_list"} as returned
#'   from \code{\link{bridge_sampler}}. Additionally, the default method assumes
#'   that all passed objects are numeric log marginal likelihoods (e.g., from
#'   \code{\link{logml}}) and will throw an error otherwise.
#' @param ... further objects of class \code{"bridge"} or \code{"bridge_list"}
#'   as returned from \code{\link{bridge_sampler}}. Or numeric values for the
#'   default method.
#' @param prior_prob numeric vector with prior model probabilities. If omitted,
#'   a uniform prior is used (i.e., all models are equally likely a priori). The
#'   default \code{NULL} corresponds to equal prior model weights.
#' @param model_names If \code{NULL} (the default) will use model names derived
#'   from deparsing the call. Otherwise will use the passed values as model
#'   names.
#'
#' @return For the default method and the method for \code{"bridge"} objects, a
#'   named numeric vector with posterior model probabilities (i.e., which sum to
#'   one).
#'
#'   For the method for \code{"bridge_list"} objects, a matrix consisting of
#'   posterior model probabilities where each row sums to one and gives the
#'   model probabilities for one set of logmls. The (named) columns correspond
#'   to the models and the number of rows is given by the \code{"bridge_list"}
#'   element with the most \code{repetitions}. Elements with fewer repetitions
#'   will be recycled (with warning).
#' @author Quentin F. Gronau and Henrik Singmann
#' @note For realistic examples, see \code{\link{bridge_sampler}} and the
#'   accompanying vignettes: \cr \code{vignette("bridgesampling_example_jags")}
#'   \cr \code{vignette("bridgesampling_example_stan")}
#' @example examples/example.post_prob.R
#' @importFrom methods is
post_prob <- function (x, ..., prior_prob = NULL, model_names = NULL) {
   UseMethod("post_prob", x)
}

#' @rdname post_prob
#' @export
post_prob.bridge <- function(x, ..., prior_prob = NULL, model_names = NULL) {
  dots <- list(...)
  mc <- match.call()
  modb <- vapply(dots, inherits, NA, what = c("bridge", "bridge_list"))
  if (is.null(model_names))
    model_names <- c(deparse(mc[["x"]]), vapply(which(modb), function(x) deparse(mc[[x+2]]), ""))
  if (sum(modb) == 0)
    stop("Only one object of class 'bridge' or 'bridge_list' passed.", call. = FALSE)
  if (sum(modb) != length(dots))
    warning("Objects not of class 'bridge' or 'bridge_list' are ignored.", call. = FALSE)

  logml <- vapply(c(list(x), dots[modb]), logml, FUN.VALUE = 0)

  .post_prob_calc(logml=logml, model_names = model_names, prior_prob=prior_prob)
}


#' @rdname post_prob
#' @export
post_prob.bridge_list <- function(x, ..., prior_prob = NULL, model_names = NULL) {
  dots <- list(...)
  mc <- match.call()
  modb <- vapply(dots, inherits, NA, what = c("bridge", "bridge_list"))
  if (is.null(model_names))
    model_names <- c(deparse(mc[["x"]]), vapply(which(modb), function(x) deparse(mc[[x+2]]), ""))
  if (sum(modb) == 0)
    stop("Only one object of class 'bridge' or 'bridge_list' passed.", call. = FALSE)
  if (sum(modb) != length(dots))
    warning("Objects not of class 'bridge' or 'bridge_list' are ignored.", call. = FALSE)

  logml <- lapply(c(list(x), dots[modb]), "[[", i = "logml")
  len <- vapply(logml, length, FUN.VALUE = 0)
  if (!all(len == max(len))) {
    warning("Not all objects provide ", max(len), " logmls. Some values are recycled.", call. = FALSE)
    logml <- lapply(logml, function(x) rep(x, length.out = max(len)))
  }
  t(apply(as.data.frame(logml), 1, .post_prob_calc,
          model_names = model_names, prior_prob=prior_prob))
}

#' @rdname post_prob
#' @export
post_prob.default <- function(x, ..., prior_prob = NULL, model_names = NULL) {
  dots <- list(...)
  mc <- match.call()
  if (is.null(model_names))
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

  if(is.null(prior_prob))
    prior_prob <- rep(1/length(logml), length(logml))

  if(!isTRUE(all.equal(sum(prior_prob), 1)))
    stop("Prior model probabilities do not sum to one.", call. = FALSE)

  if(length(logml) != length(prior_prob))
    stop("Number of objects/logml-values needs to match number of elements in prior_prob.", call. = FALSE)

  if(any(is.na(logml))) {
    post_prob <- rep(NA_real_, length(logml))
    warning("NAs in logml values. No posterior probabilities calculated.", call. = FALSE)
  } else {
    post_prob <- as.numeric(e^logml*prior_prob / sum(e^logml*prior_prob))
    if(!isTRUE(all.equal(sum(post_prob), 1)))
      warning("Posterior model probabilities do not sum to one.", call. = FALSE)
  }
  names(post_prob) <- make.unique(as.character(model_names))

  return(post_prob)

}

