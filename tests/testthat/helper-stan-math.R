# rstan and cmdstanr model methods each contain Stan Math. With GCC on Linux,
# its autodiff stack is a unique symbol that the loader binds once per
# process, so loading the second one segfaults. cmdstanr >= 1.0.0 hides the
# symbol. These guards skip a test that would load the second one.

.cmdstanr_fixed_version <- "1.0.0"

# AutodiffStackSingleton<vari_base, chainable_alloc>::instance_
.stan_math_symbol <- paste0(
  "_ZN4stan4math22AutodiffStackSingletonINS0_9vari_baseENS0_15",
  "chainable_allocEE9instance_E"
)

# cmdstanr model methods (sourceCpp_<n>) that export the autodiff stack.
# rstan's model DLLs (file<hex>) share rstan's Stan Math, so they do not count.
.cmdstanr_methods_loaded <- function() {
  any(vapply(getLoadedDLLs(), function(dll) {
    startsWith(dll[["name"]], "sourceCpp_") && !is.null(tryCatch(
      getNativeSymbolInfo(.stan_math_symbol, PACKAGE = dll),
      error = function(e) NULL
    ))
  }, logical(1)))
}

# A DLL of StanHeaders or of a package built against it (rstan, rstanarm)
.rstan_loaded <- function() {
  any(vapply(names(getLoadedDLLs()), function(name) {
    if (identical(name, "StanHeaders")) {
      return(TRUE)
    }
    linking_to <- tryCatch(
      utils::packageDescription(name, fields = "LinkingTo"),
      error = function(e) NA, warning = function(w) NA
    )
    !is.na(linking_to) && grepl("StanHeaders", linking_to, fixed = TRUE)
  }, logical(1)))
}

# Call before a test loads rstan, also before skip_if_not_installed("rstan")
skip_if_cmdstanr_methods_loaded <- function() {
  if (.cmdstanr_methods_loaded()) {
    testthat::skip("cmdstanr model methods are loaded; loading rstan would segfault")
  }
}

# Call before a test uses cmdstanr model methods
skip_if_rstan_loaded <- function() {
  fixed <- tryCatch(
    utils::packageVersion("cmdstanr") >= .cmdstanr_fixed_version,
    error = function(e) FALSE
  )
  if (!fixed && .rstan_loaded()) {
    testthat::skip("rstan is loaded; cmdstanr model methods would segfault")
  }
}
