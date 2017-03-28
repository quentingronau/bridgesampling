#' Generic function that returns log marginal likelihood from bridge objects. For objects of class \code{"bridge_list"}, which contains multiple log marginal likelihoods, \code{fun} is performed on the vector and its result returned.
#' @title Log Marginal Likelihoods from Bridge Objects
#' @param x Object of class \code{"bridge"} or \code{"bridge_list"} as returned from \code{\link{bridge_sampler}}.
#' @param fun Function which returns a scalar value and is applied to the \code{logml} vector of \code{"bridge_list"} objects. Default is \code{\link{median}}.
#' @param ... Further arguments passed to \code{fun}.
#' @return scalar numeric
#' @export
logml <- function (x, ...) {
   UseMethod("logml", x)
}


#' @rdname logml
#' @export
logml.bridge <- function (x, ...) {
   x$logml
}

#' @rdname logml
#' @export
logml.bridge_list <- function (x, fun = median, ...) {
   out <- fun(x$logml, ...)
   if (length(out) != 1) {
     warning("fun returns results of length != 1, only first used.")
     out <- out[1]
   }
   out
}

