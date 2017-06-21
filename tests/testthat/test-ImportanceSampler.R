
context("Importance Sampling")

set.seed(1234)
n <- 10
p <- 5
Niter <-  10
Group <- rep(1:p)
Weights <- rep(1, p)
x <- matrix(rnorm(n*p), n)

# Target distribution parameter
pluginTarget <- rep(0,p)
sig2Target <- .5
lbdTarget <- .37

# Proposal distribution parameter
pluginProp1 <- rep(0,p) ;pluginProp2 <- rep(1,p)
sig2Prop1 <- .5; sig2Prop2 <- 1
lbdProp1 <- .37; lbdProp2 <- .5

#
# Using non-mixture distribution
# ------------------------------
# Target distribution parameters (coeff, sig2, lbd) = (rep(0,p), .5, .37)
# Proposal distribution parameters (coeff, sig2, lbd) = (rep(1,p), 1, .5)
#
DS <- DirectSampler(X = x, pointEstimate_1 = rep(1, p), sig2_1 = 1, lbd_1 = .5,
 weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = FALSE)

test_that("Low dimensional setting", {
  expect_error(hdIS(X = x,pluginTarget = rep(0,p), sig2Target = .5, lbdTarget = .37,
                    pluginProp1 = rep(1,p),sig2Prop1 = 1,lbdProp1 = .5,proposalsample = DS, group = Group,
                    weights = Weights, log = TRUE)
  , "High dimensional setting")
})
hdIS(X = x,pluginTarget = rep(0,p), sig2Target = .5, lbdTarget = .37,
     pluginProp1 = rep(1,p),sig2Prop1 = 1,lbdProp1 = .5,proposalsample = DS, group = Group,
     weights = Weights, log = TRUE)

#
# Using mixture distribution
# ------------------------------
# Target distribution parameters (coeff, sig2, lbd) = (rep(0,p), .5, .37)
# Proposal distribution parameters
#  (coeff, sig2, lbd) = (rep(0,p), .5, .37) & (rep(1,p), 1, .5)
#
#
pluginTarget <- rep(0,p)
sig2Target <- .5
lbdTarget <- .37
pluginProp1 <- rep(0,p)
pluginProp2 <- rep(1,p)
sig2Prop1 <- .5
sig2Prop2 <- 1
lbdProp1 <- .37
lbdProp2 <- .5

DSMixture <- DirectSampler(X = x, pointEstimate_1 = pluginProp1, sig2_1 = sig2Prop1, lbd_1 = lbdProp1,
 pointEstimate_2 = pluginProp2, sig2_2 = sig2Prop2, lbd_2 = lbdProp2, weights = Weights,
 group = Group, niter = Niter, type = "coeff", parallel = TRUE)
hdIS(X = x,pluginTarget = pluginTarget, sig2Target = sig2Target, lbdTarget = lbdTarget,
     pluginProp1 = pluginProp1, sig2Prop1 = sig2Prop1, lbdProp1 = lbdProp1,
     pluginProp2 = pluginProp2, sig2Prop2 = sig2Prop2, lbdProp2 = lbdProp2,
     proposalsample = DSMixture, group = Group,
     weights = Weights, log = TRUE)
