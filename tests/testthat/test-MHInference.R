
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
  expect_error(Lasso.MHLS(X = X,Y = Y,lbd = .5)
               , NA)
  expect_error(Lasso.MHLS(X = X,Y = Y,lbd = -.5)
               , "lbd has to be positive")
  expect_error(Lasso.MHLS(X = X,Y = c(0,Y))
               , "dimension")
  expect_error(Lasso.MHLS(X = X,Y = Y,weights=c(1,1), group=rep(c(1,2),each=5))
               , NA)

  expect_error(Lasso.MHLS(X = X,Y = Y,weights=1:(p+1))
               , "length")
  expect_error(Lasso.MHLS(X = X,Y = Y,weights=1:p, group=rep(c(1,2),each=5))
               , "length")

})

set.seed(123)
n <- 20
p <- 100
X <- matrix(rnorm(n*p),n)
Y <- X %*% c(1,1,rep(0,p-2)) + rnorm(n)

test_that("High dimensional setting", {
  expect_error(Lasso.MHLS(X = X,Y = Y,lbd = .5)
               , NA)
  expect_error(Lasso.MHLS(X = X,Y = Y,lbd = -.5)
               , "lbd has to be positive")
  expect_error(Lasso.MHLS(X = X,Y = c(0,Y))
               , "dimension")
  expect_error(Lasso.MHLS(X = X,Y = Y,weights=c(1,1), group=rep(c(1,2),each=50))
               , NA)

  expect_error(Lasso.MHLS(X = X,Y = Y,weights=1:(p+1))
               , "length")
  expect_error(Lasso.MHLS(X = X,Y = Y,weights=1:p, group=rep(c(1,2),each=50))
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
LassoResult <- Lasso.MHLS(X = X,Y = Y,lbd = lbd,group=group,
weights = weights)
B0 <- LassoResult$B0
S0 <- LassoResult$S0

test_that("High dimensional setting", {
  expect_error(Postinference.MHLS(X=X, Y=Y, B0=B0, S0=S0, lbd=lbd, weights = weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
               , NA)
  expect_error(Postinference.MHLS(X=X, Y=Y, B0=B0, lbd=lbd, weights = weights, tau = rep(1, sum(B0!=0)),
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
               , "missing")
  expect_error(Postinference.MHLS(X=X, Y=Y, S0=S0, lbd=lbd, weights = weights, tau = rep(1, sum(B0!=0)),
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
               , "missing")
  expect_error(Postinference.MHLS(X=X, Y=Y, B0=B0, S0=rep(0,p), lbd=lbd, weights = weights, tau = rep(1, sum(B0!=0)),
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
               , "Invalid B0 or S0")
  expect_error(Postinference.MHLS(X=X, Y=Y, B0=rep(0,p), S0=S0, lbd=lbd, weights = weights, tau = rep(1, sum(B0!=0)),
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
               , "Invalid B0 or S0")
  expect_error(Postinference.MHLS(X=X, Y=Y, B0=B0, S0=S0, lbd=lbd, weights = -weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
               , "positive")
  expect_error(Postinference.MHLS(X=X, Y=c(Y,0), B0=B0, S0=S0, lbd=lbd, weights = weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
               , "dimension")
  expect_error(Postinference.MHLS(X=X[-1,], Y=Y, B0=B0, S0=S0, lbd=lbd, weights = weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
               , "dimension")
  expect_error(Postinference.MHLS(X=X[,-1], Y=Y, B0=B0, S0=S0, lbd=lbd, weights = weights,
                                  sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
               , "length")

  expect_warning(Postinference.MHLS(X=X, Y=Y, B0=B0, S0=S0, lbd=lbd, weights = weights,
                                    sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE, ncores = 10)
                 , "ncores is larger")
  expect_warning(Postinference.MHLS(X=X, Y=Y, B0=B0, S0=S0, lbd=lbd, weights = weights,
                                    sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE, ncores = 1)
                 , "needs to be greater than 1")
})


