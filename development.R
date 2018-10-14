require(devtools)
require(testthat)
options(error = NULL)

load_all()
devtools::test()

roxygen2::roxygenize()

build_vignettes()
devtools::build()

Sys.setenv(`_R_CHECK_FORCE_SUGGESTS_` = "false")
devtools::check()

### check reverse dependencies:

#revdep()
devtools::revdep_check(libpath = "../revdep", check_dir = "../revdep_checks")
#devtools::install.packages("brms", lib = "../revdep")
devtools::revdep_check_resume()
devtools::revdep_check_save_summary()
devtools::revdep_check_print_problems()
devtools::revdep_maintainers()
