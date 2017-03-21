require(devtools)
require(testthat)
options(error = NULL)

load_all()
test()
test_file(path = "test-bridge_sampler_normal_Rcpp_parallel.R")

roxygen2::roxygenize()
