#' Methods for bridge and bridge_list objects
#'
#' Methods defined for objects returned from the generic \code{\link{bridge_sampler}} function.
#'
#' @param object,x object of class \code{bridge} or \code{bridge_list} as returned from \code{\link{bridge_sampler}}.
#' @param na.rm logical. Should NA estimates in \code{bridge_list} objects be removed? Passed to \code{\link{error_measures}}.
#' @param ... further arguments, currently ignored.
#'
#' @return
#' The \code{summary} methods return a \code{data.frame} which contains the log marginal likelihood plus the result returned from invoking \code{\link{error_measures}}.
#'
#' The \code{print} methods simply print and return nothing.
#'
#'
#' @name bridge-methods
NULL


# summary methods

#' @rdname bridge-methods
#' @method summary bridge
#' @export
summary.bridge <- function(object, na.rm = TRUE, ...) {

  if( ! (object$method %in% c("normal", "warp3"))) {
    stop('object$method needs to be either "normal" or "warp3".', call. = FALSE)
  }

  if (object$method == "normal") {

    em <- error_measures(object)
    out <- data.frame("Logml_Estimate" = object$logml,
                      "Relative_Mean_Squared_Error" = em$re2,
                      "Coefficient_of_Variation" = em$cv,
                      "Percentage_Error" = em$percentage,
                      "Method" = object$method,
                      "Repetitions" = 1,
                      stringsAsFactors = FALSE)

  } else if (object$method == "warp3") {

    out <- data.frame("Logml_Estimate" = object$logml,
                      "Method" = object$method,
                      "Repetitions" = 1)

  }

  class(out) <- c("summary.bridge", "data.frame")
  return(out)

}

#' @rdname bridge-methods
#' @method summary bridge_list
#' @export
summary.bridge_list <- function(object, na.rm = TRUE, ...) {

  if( ! (object$method %in% c("normal", "warp3"))) {
    stop('object$method needs to be either "normal" or "warp3".', call. = FALSE)
  }

  em <- error_measures(object, na.rm = na.rm)
  out <- data.frame("Logml_Estimate" = median(object$logml, na.rm = na.rm),
                    "Min" = em$min,
                    "Max" = em$max,
                    "Interquartile_Range" = em$IQR,
                    "Method" = object$method,
                    "Repetitions" = object$repetitions)

  class(out) <- c("summary.bridge_list", "data.frame")
  return(out)

}

# print summary methods

#' @rdname bridge-methods
#' @method print summary.bridge
#' @export
print.summary.bridge <- function(x, ...) {

  if (x[["Method"]] == "normal") {

    cat('\nBridge sampling log marginal likelihood estimate \n(method = "',
        as.character(x[["Method"]]),
        '", repetitions = ', x[["Repetitions"]], '):\n\n ',
          x[["Logml_Estimate"]],
        '\n\nError Measures:\n\n Relative Mean-Squared Error: ',
          x[["Relative_Mean_Squared_Error"]],
        '\n Coefficient of Variation: ', x[["Coefficient_of_Variation"]],
        '\n Percentage Error: ', x[["Percentage_Error"]],
        '\n\nNote:\nAll error measures are approximate.\n\n', sep = "")

  } else if (x[["Method"]] == "warp3") {

    cat('\nBridge sampling log marginal likelihood estimate \n(method = "',
        as.character(x[["Method"]]),
        '", repetitions = ', x[["Repetitions"]], '):\n\n ',
        x[["Logml_Estimate"]],
        '\n\nNote:\nNo error measures are available for method = "warp3"',
        '\nwith repetitions = 1.',
        '\nWe recommend to run the warp3 procedure multiple times to',
        '\nassess the uncertainty of the estimate.\n\n', sep = "")

  }

}

#' @rdname bridge-methods
#' @method print summary.bridge_list
#' @export
print.summary.bridge_list <- function(x, ...) {
  cat('\nBridge sampling log marginal likelihood estimate \n(method = "',
      as.character(x[["Method"]]), '", repetitions = ', x[["Repetitions"]],
      '):\n\n ', x[["Logml_Estimate"]],
      '\n\nError Measures:\n\n Min: ', x[["Min"]],
      '\n Max: ', x[["Max"]],
      '\n Interquartile Range: ', x[["Interquartile_Range"]],
      '\n\nNote:\nAll error measures are based on ', x[["Repetitions"]],
      ' estimates.\n\n', sep = "")

}

# print methods

#' @rdname bridge-methods
#' @method print bridge
#' @export
print.bridge <- function(x, ...) {

  cat("Bridge sampling estimate of the log marginal likelihood: ",
      round(x$logml, 5), "\nEstimate obtained in ", x$niter,
      " iteration(s) via method \"", x$method, "\".\n", sep = "")
}

#' @rdname bridge-methods
#' @method print bridge_list
#' @export
print.bridge_list <- function(x, na.rm = TRUE, ...) {

  cat("Median of ", x$repetitions,  " bridge sampling estimates\nof the log marginal likelihood: ",
      round(median(x$logml, na.rm = na.rm), 5), "\nRange of estimates: ", round(range(x$logml, na.rm = na.rm)[1], 5), " to ",
      round(range(x$logml, na.rm = na.rm)[2], 5),
      "\nInterquartile range: ", round(stats::IQR(x$logml, na.rm = na.rm), 5), "\nMethod: ", x$method,  "\n", sep = "")
  if (any(is.na(x$logml))) warning(sum(is.na(x$logml))," bridge sampling estimate(s) are NAs.", call. = FALSE)
}


