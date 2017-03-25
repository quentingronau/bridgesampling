require(devtools)
require(testthat)
options(error = NULL)

load_all()
test()

roxygen2::roxygenize()

# use_vignette("hierarchical_normal_example")
# build_vignettes()
