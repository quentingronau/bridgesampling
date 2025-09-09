<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/bridgesampling)](https://CRAN.R-project.org/package=bridgesampling)
[![R-CMD-check](https://github.com/quentingronau/bridgesampling/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/quentingronau/bridgesampling/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

bridgesampling: Bridge Sampling for Marginal Likelihoods and Bayes Factors
====

`bridgesampling` is an R package for conducting Bayesian model comparisons using bridge sampling (Meng & Wong, 1996).
Specifically, it allows one to compute marginal likelihoods, Bayes factors, and posterior model probabilities.

Meng, X.-L., & Wong, W. H. (1996). Simulating ratios of normalizing constants via a simple identity: A theoretical exploration. *Statistica Sinica*, 6, 831-860.

For additional information, see the vignettes with a hierarchical normal example implemented in Jags ([link](https://htmlpreview.github.io/?https://github.com/quentingronau/bridgesampling/blob/master/inst/doc/bridgesampling_example_jags.html)) and in Stan ([link](https://htmlpreview.github.io/?https://github.com/quentingronau/bridgesampling/blob/master/inst/doc/bridgesampling_example_stan.html)), and a Bayesian t-test example implemented in Stan ([link](https://htmlpreview.github.io/?https://github.com/quentingronau/bridgesampling/blob/master/inst/doc/bridgesampling_stan_ttest.html)).


## Installation

- `bridgesampling` is [available from CRAN](https://cran.r-project.org/package=bridgesampling) so the current stable version can be installed directly via: `install.packages("bridgesampling")`

- To install the latest development version you will need the [`devtools`](https://github.com/r-lib/devtools) package: 
  `devtools::install_github("quentingronau/bridgesampling@master")`
  
- For building, use `--compact-vignettes="gs+qpdf"`

----
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
