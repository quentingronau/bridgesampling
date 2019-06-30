## ------------------------------------------------------------------------
library(bridgesampling)

### generate data ###
set.seed(12345)

mu <- 0
tau2 <- 0.5
sigma2 <- 1

n <- 20
theta <- rnorm(n, mu, sqrt(tau2))
y <- rnorm(n, theta, sqrt(sigma2))
  

## ----eval=FALSE----------------------------------------------------------
#  ### set prior parameters ###
#  mu0 <- 0
#  tau20 <- 1
#  alpha <- 1
#  beta <- 1

## ---- eval=FALSE---------------------------------------------------------
#  library("nimble")
#  
#  # models
#  codeH0 <- nimbleCode({
#    invTau2 ~ dgamma(1, 1)
#    tau2 <- 1/invTau2
#    for (i in 1:20) {
#      theta[i] ~ dnorm(0, sd = sqrt(tau2))
#      y[i] ~ dnorm(theta[i], sd = 1)
#    }
#  })
#  codeH1 <- nimbleCode({
#    mu ~ dnorm(0, sd = 1)
#    invTau2 ~ dgamma(1, 1)
#    tau2 <- 1/invTau2
#    for (i in 1:20) {
#      theta[i] ~ dnorm(mu, sd = sqrt(tau2))
#      y[i] ~ dnorm(theta[i], sd = 1)
#    }
#  })
#  
#  ## steps for H0:
#  modelH0 <- nimbleModel(codeH0)
#  modelH0$setData(y = y) # set data
#  cmodelH0 <- compileNimble(modelH0) # make compiled version from generated C++
#  
#  ## steps for H1:
#  modelH1 <- nimbleModel(codeH1)
#  modelH1$setData(y = y) # set data
#  cmodelH1 <- compileNimble(modelH1) # make compiled version from generated C++
#  

## ---- eval=FALSE---------------------------------------------------------
#  
#  # build MCMC functions, skipping customization of the configuration.
#  mcmcH0 <- buildMCMC(modelH0,
#                      monitors = modelH0$getNodeNames(stochOnly = TRUE,
#                                                      includeData = FALSE))
#  mcmcH1 <- buildMCMC(modelH1,
#                      monitors = modelH1$getNodeNames(stochOnly = TRUE,
#                                                      includeData = FALSE))
#  # compile the MCMC function via generated C++
#  cmcmcH0 <- compileNimble(mcmcH0, project = modelH0)
#  cmcmcH1 <- compileNimble(mcmcH1, project = modelH1)
#  
#  # run the MCMC.  This is a wrapper for cmcmc$run() and extraction of samples.
#  # the object samplesH1 is actually not needed as the samples are also in cmcmcH1
#  samplesH0 <- runMCMC(cmcmcH0, niter = 1e5, nburnin = 1000, nchains = 2,
#                       progressBar = FALSE)
#  samplesH1 <- runMCMC(cmcmcH1, niter = 1e5, nburnin = 1000, nchains = 2,
#                       progressBar = FALSE)

## ---- echo=FALSE---------------------------------------------------------
load(system.file("extdata/", "vignette_example_nimble.RData",
                     package = "bridgesampling"))

## ----eval=FALSE----------------------------------------------------------
#  # compute log marginal likelihood via bridge sampling for H0
#  H0.bridge <- bridge_sampler(cmcmcH0, silent = TRUE)
#  
#  # compute log marginal likelihood via bridge sampling for H1
#  H1.bridge <- bridge_sampler(cmcmcH1, silent = TRUE)

## ------------------------------------------------------------------------
print(H0.bridge)
print(H1.bridge)

## ----eval=FALSE----------------------------------------------------------
#  # compute percentage errors
#  H0.error <- error_measures(H0.bridge)$percentage
#  H1.error <- error_measures(H1.bridge)$percentage

## ------------------------------------------------------------------------
print(H0.error)
print(H1.error)

## ------------------------------------------------------------------------
# compute Bayes factor
BF01 <- bf(H0.bridge, H1.bridge)
print(BF01)

## ------------------------------------------------------------------------
# compute posterior model probabilities (assuming equal prior model probabilities)
post1 <- post_prob(H0.bridge, H1.bridge)
print(post1)

## ------------------------------------------------------------------------
# compute posterior model probabilities (using user-specified prior model probabilities)
post2 <- post_prob(H0.bridge, H1.bridge, prior_prob = c(.6, .4))
print(post2)

