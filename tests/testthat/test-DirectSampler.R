
context("Direct Sampling")

set.seed(1234)
n <- 10
p <- 30
Niter <-  10
Group <- rep(1:(p/10), each = 10)
Weights <- rep(1, p/10)
x <- matrix(rnorm(n*p), n)
Y <- x%*%rep(1,p) + rnorm(n)

test_that("pointEstimate length under two types, \"mu\" and \"coeff\"", {
 expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p-1), sig2_1 = 1, lbd_1 = .5,
                            weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = FALSE),
              "pointEstimate must have a same length")
 expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = FALSE)
  , NA)
 expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, n-1), sig2_1 = 1, lbd_1 = .5,
                            weights = Weights, group = Group, niter = Niter, type = "mu", parallel = FALSE),
              "pointEstimate must have a same length")
 expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, n), sig2_1 = 1, lbd_1 = .5,
                            weights = Weights, group = Group, niter = Niter, type = "mu", parallel = FALSE)
              , NA)
 expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, n), sig2_1 = 1, lbd_1 = .5,
                            pointEstimate_2 = rep(0, n), sig2_2 = 1, lbd_2 = .5,
                            weights = Weights, group = Group, niter = Niter, type = "mu", parallel = FALSE)
              , NA)
})

test_that("missing parameter", {
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1,
                             weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = FALSE)
               , "provide all the parameters")
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, n), sig2_1 = 1, lbd_1 = .5,
                             pointEstimate_2 = rep(0, n), sig2_2 = 1,
                             weights = Weights, group = Group, niter = Niter, type = "mu", parallel = FALSE)
               , "provide all the parameters")
})

test_that("Improper value of parameter, sig2/lbd/weights", {
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = -1, lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, parallel = FALSE),
               "sig2 should be positive")
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = -1,
                             weights = Weights, group = Group, niter = Niter, parallel = FALSE),
               "lbd should be non-negative")
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = c(1,2), lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, parallel = FALSE),
               "sig2/lbd should be a scalar")
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = c(-1, rep(1, p / 10 - 1)), group = Group, niter = Niter, parallel = FALSE),
               "weights should be positive")
})

test_that("group argument",{
  Group <- rep(1:(p/10), each = 10); Weights <- rep(1, p/10)
  Group1 <- rep(1,p); Weights1 <- 1   # valid
  Group2 <- rep(1,p-1); Weights2 <- 1 # length(group) != p
  Group3 <- rep(c(1,2,4),each=10); Weights3 <- rep(1, 3) # group index are not consecutive integers
  Group4 <- rep(-1,p); Weights4 <- 1  # using negative integer
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = FALSE)
               , NA)
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights1, group = Group1, niter = Niter, type = "coeff", parallel = FALSE)
  , NA)
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights2, group = Group2, niter = Niter, type = "coeff", parallel = FALSE)
  , "group must have a same length with the number of X columns")
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights3, group = Group3, niter = Niter, type = "coeff", parallel = FALSE)
  , "group index has to be a consecutive integer starting from 1"
  )
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights4, group = Group4, niter = Niter, type = "coeff", parallel = FALSE)
  , "group index has to be a consecutive integer starting from 1"
  )
})

test_that("method is either nonparametric or normal", {
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, type = "coeff", method = "nonparametric", parallel = FALSE)
               , "Y is needed")
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, Y=Y, type = "coeff", method = "nonparametric", parallel = FALSE)
               , NA)
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, Y=Y, type = "coeff", method = "random", parallel = FALSE)
               , "either normal or nonparametric")
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, Y=Y, type = "coeff", method = "normal", parallel = FALSE)
               , NA)
  expect_error(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, type = "coeff", method = "normal", parallel = FALSE)
               , NA)
})

test_that("parallel", {
  expect_warning(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                             weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = TRUE,
                             ncores = 1)
               , "Set ncores to 2")
  # expect_warning(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
  #                              weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = TRUE,
  #                              ncores = 100000)
  #                , "maximum possible value")
  expect_warning(DirectSampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
                               weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = TRUE,
                               ncores = 2)
                , NA)
})


