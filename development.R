require(devtools)
require(testthat)
load_all()
test()

roxygen2::roxygenize()
