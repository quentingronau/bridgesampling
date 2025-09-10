# bridgesampling 1.2-0 (2025-09-10)

* Added CmdStanR method and corresponding tests (thanks to @GiorgioMB and @avehtari #44).
* Added Monte Carlo Standard Error (MCSE) to bridgesampling, see: https://arxiv.org/abs/2508.14487 (thanks to @GiorgioMB and @avehtari #43).
* Fixed bug in simplex with small dimensionality (thanks to @FBartos #31).
* Added a `NEWS.md` file to track changes to the package.

# bridgesampling 1.1-5 (2023-06-01)

* Deactivated stanreg tests to avoid CRAN check issues.

# bridgesampling 1.1-0 (2021-03-01)

* Fixed subscript out of bounds error, see: https://github.com/quentingronau/bridgesampling/issues/26
* Deactivated stan tests on Windows to avoid CRAN check issues.

# bridgesampling 1.0-0 (2020-02-01)

* Included citation file and references to JSS article

# bridgesampling 0.8-0 (2019-12-01)

* Disabled use of mvnfast and revetred back to mvtnorn. see also: https://github.com/quentingronau/bridgesampling/issues/20
* Version 0.7-x introduced a bug that prevented a rerunning of the iterative scheme based on harmonic mean in case maxit was reached. This bug should now be removed. See: https://github.com/quentingronau/bridgesampling/issues/18
* For older news see NEWS.old file on GitHub: https://github.com/quentingronau/bridgesampling/blob/master/NEWS.old

