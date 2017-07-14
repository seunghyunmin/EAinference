#' @title lasso / group lasso estimator
#'
#' @description provides lasso / group lasso solutio;
#' coefficient-estimate and subgradient.
#'
#' @param X n x p matrix of predictors.
#' @param Y n x 1 vector of response.
#' @param lbd penalty term of lasso. See the loss function given below for more
#'  details.
#' @param weights Weight term for each group. Default is
#' \code{rep(1, max(group))}.
#' @param group p x 1 vector of consecutive integers. The number of groups
#' should be same as max(group).
#' @details
#' Using gglasso package, provide lasso / group lasso solution along with
#' subgradient. The loss function for group lasso is
#' \deqn{L(\beta) = ||y-X\beta||^2 / (2n) + \lambda \sum_j ||\beta_(j)||,}
#' where (j) is the index set of j-th group. If the size of the group is 1 for
#' every group, it becomes lasso loss function.
#' @return \item{B0}{a vector of coefficient estimator.}
#' @return \item{S0}{a vector of subgradient.}
#' @return \item{lbd, weights, grouop}{same as input arguments.}
#' @examples
#' set.seed(123)
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n*p),n)
#' Y <- X %*% c(1,1,rep(0,p-2)) + rnorm(n)
#' #
#' # lasso
#' #
#' Lasso.MHLS(X = X,Y = Y,lbd = .5)
#' #
#' # group lasso
#' #
#' Lasso.MHLS(X = X,Y = Y,lbd = .5,weights = rep(1,2),group=rep(1:2,each=5))
#' @export
Lasso.MHLS <- function(X, Y, lbd=.37, weights=rep(1,max(group)),
  group=1:ncol(X))
{
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)
  Y <- matrix(Y, , 1)

  if (lbd <= 0) {
    stop("lbd has to be positive.")
  }
  if (length(group) != p) {
    stop("length(group) has to be the same with ncol(X)")
  }
  if (length(weights) != length(unique(group))) {
    stop("length(weights) has to be the same with the number of groups")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (any(!1:max(group) %in% group)) {
    stop("group index has to be a consecutive integer starting from 1.")
  }

  if (length(Y) != n) {
    stop("dimension of X and Y are not conformable.")
  }



  IndWeights <- rep(weights,table(group))
  # scale X with weights
  X.tilde   <- scale(X,FALSE,scale=IndWeights)
  # compute group lasso estimator B0 and S0
  B0 <- coef(gglasso(X.tilde, Y, pf = rep(1,max(group)), group = group,
                     loss="ls", intercept=F, lambda=lbd))[-1] / IndWeights
  S0 <- (t(X.tilde) %*% Y - t(X.tilde) %*% X.tilde %*%
           (B0 * IndWeights)) / n / lbd
  #A <- which(B0!=0)
  return(list(B0=B0, S0=c(S0), lbd=lbd, weights=weights, group=group))
}

#' @title Post-inference for lasso estimator
#'
#' @description Provides Confidence intervals for the set of active coefficients
#' from lasso estimator using Metropolis-Hastings sampler.
#'
#' @param X predictors matrix.
#' @param Y response vector.
#' @param lbd penalty term of lasso. See the loss function given below for
#' more details.
#' @param weights Weight term for each group. Default is
#' \code{rep(1, max(group))}.
#' @param tau numeric vector. Standard deviaion of proposal distribution
#'  for each active beta.
#' Adjust the value to get relevant level of acceptance rate.
#' @param sig2.hat variance of error term.
#' @param alpha confidence level for confidence interval.
#' @param nChain The number of chains each from different pluginbeta.
#' @param niterPerChain The number of iteration per chain.
#' @param parallel Whether to parallelize the code. Default is \code{FALSE}.
#' @param ncores The number of cores to use for the parallelization. If missing,
#'  it uses maximum number of cores.
#' @param printSamples Boolean value. If true, print Metropolis-Hastings samples.
#' @param ... auxiliary \code{link{MHLS}} arguments.
#' @details
#' We provide Post-selection inference for lasso estimator.
#' Using Metropolis Hastings Sampler with multiple chanin, \code{(1-alpha)}
#' confidence interval for each active coefficients is generated.
#' Set \code{printSamples=TRUE} if one wants to check the samples.
#' Check the acceptance rate and adjust tau accordingly. Desirable level of
#' acceptance rate for beta is \code{0.30 +- 0.15}. We recommend to set
#' nChain >= 10 and niterPerChain >=500.
#' @return \item{MHsamples}{a list of a class MHLS.}
#' @return \item{confidenceInterval}{(1-alpha) confidence interval
#' of each active coefficient.}
#' @examples
#' set.seed(123)
#' n <- 5
#' p <- 10
#' X <- matrix(rnorm(n*p),n)
#' Y <- X %*% rep(1,p) + rnorm(n)
#' sig2 <- 1
#' lbd <- .37
#' weights <- rep(1,p)
#' group <- 1:p
#' Postinference.MHLS(X, Y, lbd, weights, tau=rep(1, sum(B0!=0)),
#' sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20, parallel=TRUE)
#' Postinference.MHLS(X, Y, lbd, weights, tau=rep(1, sum(B0!=0)),
#' sig2.hat=1, alpha=.05, nChain=3, niterPerChain=20,
#' parallel=TRUE, printSamples=TRUE)
#' @export
Postinference.MHLS <- function(X, Y, lbd, weights = rep(1, length(B0)),
  tau = rep(1, sum(B0!=0)), sig2.hat, alpha = .05, nChain = 10,
  niterPerChain = 500, parallel = FALSE, ncores = 2L, printSamples=FALSE, ...)
{
  # nChain : the number of MH chains
  # niterPerChain : the number of iteration for each chain
  # B0, S0 : The lasso estimator
  # tau : same as in MHLS function

  LassoEst <- Lasso.MHLS(X, Y, lbd=lbd, weights=weights)
  B0 <- LassoEst$B0
  S0 <- LassoEst$S0
  A <- which(B0!=0)

  if (length(A)==0) {
    stop("Given lbd, active set is empty.")
  }

  if(.Platform$OS.type == "windows" && parallel == TRUE){
    n.cores <- 1L
    parallel <- FALSE
    warning("Under Windows platform, parallel computing cannot be executed.")
  }

  Y <- matrix(Y, , 1)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  nChain <- as.integer(nChain)
  niterPerChain <- as.integer(niterPerChain)

  #--------------------
  # Error Handling
  #--------------------
  if (nrow(X) != nrow(Y)) {
    stop("The dimension of X and Y are not conformable.")
  }
  if (sig2.hat <=0 || lbd <= 0) {
    stop("sig2.hat and/or lbd have to be positive.")
  }
  if (length(weights) != ncol(X)) {
    stop("length(weights) has to be the same with col(X).")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (alpha <=0 || alpha >=1) {
    stop("alpha needs to be between 0 and 1.")
  }
  if (parallel && ncores == 1) {
    ncores <- 2
    warning("If parallel=TRUE, ncores needs to be greater than 1. Automatically
Set ncores to 2.")
  }
  if (parallel && (ncores > parallel::detectCores())) {
    ncores <- parallel::detectCores()
    warning("ncores is larger than the maximum number of available processes.
Set it to the maximum possible value.")
  }
  if (any(c(nChain,niterPerChain) <= 0)) {
    stop("nChain & niterPerChain have to be a positive integer.")
  }
  # if (all.equal(coef(gglasso(X, Y, pf = weights, group = 1:p, loss="ls",
  #                             intercept=F, lambda=lbd))[-1],B0) != TRUE ||
  #     all.equal(c(((t(X)/weights)%*%Y - (t(X) /weights) %*% X %*% B0) / n / lbd)
  #                , S0) != TRUE) {
  #   stop("Invalid B0 or S0, use Lasso.MHLS to get a valid lasso solution.")
  # }
  # Draw samples of pluginbeta from the 95% confidence
  #  region boundary of restricted lse.
  # If nChain ==1, we just use restricted lse.
  if (missing(sig2.hat)) {
    if (length(A) >= nrow(X)) {
      stop("If size of active set matches with nrow(X), sig2.hat needs to be provided.")
    }
    sig2.hat <- summary((lm(Y~X[,A]+0)))$sigma^2
  }

  Pluginbeta.seq <- Pluginbeta.MHLS(X,Y,A,nChain,sqrt(sig2.hat))

  FF <- function(x) {
    MHLS(X = X, pointEstimate = Pluginbeta.seq[x,], sig2 = sig2.hat, lbd = lbd,
         weights = weights, niter=niterPerChain,
         burnin = 0, B0 = B0, S0 = S0, tau = tau, type = "coeff", verbose=FALSE, ...)
  }

  if (parallel) {
    TEMP <- parallel::mclapply(1:nChain,FF, mc.cores = ncores)
  } else {
    TEMP <- lapply(1:nChain,FF)
  }
  names(TEMP) <- paste("Chain",1:nChain,sep="")

  MCSAMPLE <- TEMP[[1]]
  if (nChain > 1) {
    for (i in 2:nChain) {
      MCSAMPLE$beta <- rbind(MCSAMPLE$beta, TEMP[[i]]$beta)
      MCSAMPLE$subgrad <- rbind(MCSAMPLE$subgrad, TEMP[[i]]$subgrad)
      MCSAMPLE$acceptHistory <- MCSAMPLE$acceptHistory + TEMP[[i]]$acceptHistory
      MCSAMPLE$niteration <- MCSAMPLE$niteration + TEMP[[i]]$niteration
      MCSAMPLE$burnin <- MCSAMPLE$burnin + TEMP[[i]]$burnin
    }
  }

  # Using MH samples, refit the coeff.
  RefitBeta <- Refit.MHLS(X,weights,lbd,MCSAMPLE)
  if (printSamples) {
    return(list(MHsamples=TEMP,pluginbeta=Pluginbeta.seq,
                          confidenceInterval=CI.MHLS(betaRefit = RefitBeta,
                          pluginbeta = Pluginbeta.seq, alpha=alpha)))
  } else {
    return(CI.MHLS(betaRefit = RefitBeta, pluginbeta = Pluginbeta.seq,
                   alpha=alpha))
  }
}

# Refit the beta estimator to remove the bias
Refit.MHLS <- function(X,weights,lbd,MHsample) {
  # W : diag(weights)
  # MHsample : MH samples from MHLS function
  n <- nrow(X)
  # Active set
  A <- which(MHsample$beta[1,]!=0)
  # Recover y using KKT condition
  Y.MH <- solve(X%*%t(X))%*%X %*% (crossprod(X) %*% t(MHsample$beta) / n +
                                     lbd * weights * t(MHsample$subgrad)) * n
  # Refit y to the restricted set of X
  beta.refit <- solve(t(X[,A])%*%X[,A])%*%t(X[,A]) %*% Y.MH
  return(beta.refit)
}

# Generate 1-alpha Confidence Interval based on the deviation
CI.MHLS <- function(betaRefit, pluginbeta, alpha=.05) {
  # pluginbeta : a nPlugin x |A| matrix of pluginbeta,
  #  note that the size is not nPlugin x p.
  # beta.refit : refitted beta via Refit.MHLS, a niter x |A| matrix.
  # alpha : significant level.
  A <- which(pluginbeta[1,]!=0)
  Quantile <- apply(betaRefit - pluginbeta[1,A], 1, function(x)
    {quantile(x,prob=c(alpha/2, 1 - alpha/2))})
  Result <- rbind(LowerBound = -Quantile[2,] + pluginbeta[1,A] ,
                  UpperBound = -Quantile[1,] + pluginbeta[1,A])
  colnames(Result) <- paste("beta", A, sep="")
  return(Result)
}

# Generate pluginbeta's from 95% confidence region
Pluginbeta.MHLS <- function(X,Y,A,nPlugin,sigma.hat) {
  # nPlugin : number of pluginbeta's want to generate
  # sigma.hat : estimator of sigma , \epsilon ~ N(0, sigma^2)
  #             If missing, use default way to generate it
  if (length(A)==0) stop("The lasso solution has an empty active set.")

  if (missing(sigma.hat)) {
    sigma.hat <- summary((lm(Y~X[,A]+0)))$sigma
  }
  beta.refit <- coef(lm(Y~X[,A]+0))

  if (nPlugin == 1) {
    return(matrix(beta.refit,1))
  } else {
    xy <- matrix(rnorm(length(A)*(nPlugin-1)), nPlugin-1)
    lambda <- 1 / sqrt(rowSums(xy^2))
    xy <- xy * lambda * sqrt(qchisq(0.95, df=length(A)))
    coeff.seq <- matrix(0,nPlugin,ncol(X))
    coeff.seq[,A]  <- rbind(beta.refit,t(t(xy%*%chol(solve(crossprod(X[,A])))) *
                                           sigma.hat + beta.refit))
    return(coeff.seq=coeff.seq)
  }
}


