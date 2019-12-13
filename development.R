require(devtools)
require(testthat)
options(error = NULL)

load_all()
devtools::test()

devtools::document()

build_vignettes()
devtools::build(args = '--compact-vignettes=gs+qpdf')

Sys.setenv(`_R_CHECK_FORCE_SUGGESTS_` = "false")
devtools::check()

### check reverse dependencies:

usethis::use_revdep()
revdepcheck::revdep_check(num_workers = 4)
