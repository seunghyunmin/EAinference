
context("Importance Sampling")

set.seed(1234)
n <- 10
p <- 5
Niter <-  10
Group <- rep(1:p)
Weights <- rep(1, p)
x <- matrix(rnorm(n*p), n)

# Target distribution parameter
pETarget <- rep(0,p)
sig2Target <- .5
lbdTarget <- .37

# Proposal distribution parameter
pEProp1 <- rep(0, p) ;pEProp2 <- rep(1, p)
sig2Prop1 <- .5; sig2Prop2 <- 1
lbdProp1 <- .37; lbdProp2 <- .5

DS <- DirectSampler(X = x, pointEstimate_1 = rep(1, p), sig2_1 = 1, lbd_1 = .5,
 weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = FALSE)

test_that("Low dimensional setting", {
  expect_error(hdIS(DirectSample = DS,pETarget = rep(0,p), sig2Target = .5, lbdTarget = .37,
                    log = TRUE)
  , "High dimensional setting")
})

class(DS) <- "list"

test_that("Wrong class", {
  expect_error(hdIS(DirectSample = DS,pETarget = rep(0,p), sig2Target = .5, lbdTarget = .37,
                    log = TRUE)
               , "Use EAlasso::DirectSampler")
})


set.seed(1234)
n <- 10
p <- 30
Niter <-  10
Group <- rep(1:(p/10), each = 10)
Weights <- rep(1, p/10)
x <- matrix(rnorm(n*p), n)

DSMixture <- DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
              pointEstimate_2 = rep(1, p), sig2_2 = 2, lbd_2 = .3, weights = Weights,
              group = Group, niter = Niter, parallel = FALSE)
test_that("Mixture", {
  expect_error(hdIS(DirectSample = DSMixture, pETarget = rep(0,p),
    sig2Target = .5, lbdTarget = .37, log = TRUE), NA)
})
