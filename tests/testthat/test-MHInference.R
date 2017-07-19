
context("MHInference")

#-------------------
# Lasso.MHLS
#-------------------
set.seed(123)
n <- 50
p <- 10
X <- matrix(rnorm(n*p),n)
Y <- X %*% c(1,1,rep(0,p-2)) + rnorm(n)

test_that("Low dimensional setting", {
  expect_error(Lasso.MHLS(X = X,Y = Y, type="lasso", lbd = .5)
               , NA)
  expect_error(Lasso.MHLS(X = X,Y = Y, type="lasso", lbd = -.5)
               , "lbd has to be positive")
  expect_error(Lasso.MHLS(X = X,Y = c(0,Y))
               , "dimension")
  expect_error(Lasso.MHLS(X = X,Y = Y, type="grlasso", weights=c(1,1), group=rep(c(1,2),each=5))
               , NA)

  expect_error(Lasso.MHLS(X = X,Y = Y, type="lasso", weights=1:(p+1))
               , "length")
  expect_error(Lasso.MHLS(X = X,Y = Y, type="grlasso", weights=1:p, group=rep(c(1,2),each=5))
               , "length")

})

set.seed(123)
n <- 20
p <- 100
X <- matrix(rnorm(n*p),n)
Y <- X %*% c(1,1,rep(0,p-2)) + rnorm(n)

test_that("High dimensional setting", {
  expect_error(Lasso.MHLS(X = X,Y = Y, type="lasso",lbd = .5)
               , NA)
  expect_error(Lasso.MHLS(X = X,Y = Y, type="lasso",lbd = -.5)
               , "lbd has to be positive")
  expect_error(Lasso.MHLS(X = X,Y = c(0,Y), type="lasso")
               , "dimension")
  expect_error(Lasso.MHLS(X = X,Y = Y, type="grlasso",weights=c(1,1), group=rep(c(1,2),each=50))
               , NA)

  expect_error(Lasso.MHLS(X = X,Y = Y, type="lasso",weights=1:(p+1))
               , "length")
  expect_error(Lasso.MHLS(X = X,Y = Y, type="grlasso",weights=1:p, group=rep(c(1,2),each=50))
               , "length")

})


#-------------------
# Postinference.MHLS
#-------------------
set.seed(123)
n <- 5
p <- 10
X <- matrix(rnorm(n*p),n)
Y <- X %*% rep(1,p) + rnorm(n)
sig2 <- 1
lbd <- .37
weights <- rep(1,p)
group <- 1:p

test_that("High dimensional setting", {
  expect_error(Postinference.MHLS(X=X, Y=Y, lbd=lbd, weights = weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel = FALSE)
               , NA)
  expect_error(Postinference.MHLS(X=X, Y=Y, lbd=lbd, weights = -weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel = FALSE)
               , "positive")
  expect_error(Postinference.MHLS(X=X, Y=c(Y,0), lbd=lbd, weights = weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel = FALSE)
               , "dimension")
  expect_error(Postinference.MHLS(X=X[-1,], Y=Y, lbd=lbd, weights = weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel = FALSE)
               , "dimension")
  expect_error(Postinference.MHLS(X=X[,-1], Y=Y, lbd=lbd, weights = weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel = FALSE)
               , "length")

  # expect_warning(Postinference.MHLS(X=X, Y=Y, B0=B0, S0=S0, lbd=lbd, weights = weights,
  #                                   sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel = FALSE, ncores = 10000)
  #                , "ncores is larger")
  if(.Platform$OS.type != "windows"){
    expect_warning(Postinference.MHLS(X=X, Y=Y, lbd=lbd, weights = weights,
                                      sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel = TRUE, ncores = 1)
                   , "needs to be greater than 1")
  } else {
    expect_warning(Postinference.MHLS(X=X, Y=Y, lbd=lbd, weights = weights,
                                      sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel = TRUE)
                   , "Under Windows platform")
  }
})



