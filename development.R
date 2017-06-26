require(devtools)
require(testthat)
options(error = NULL)

load_all()
test()

roxygen2::roxygenize()

build_vignettes()
