.onLoad <- function(libname, pkgname) {
  if (requireNamespace("posterior", quietly = TRUE)) {
    options(bridgesampling.ess_function = "posterior")
  } else {
    options(bridgesampling.ess_function = "coda")
    packageStartupMessage(
      "Package 'posterior' not found; defaulting to ESS via coda::effectiveSize(). ",
      "Install 'posterior' to use posterior-based ESS."
    )
  }
}
