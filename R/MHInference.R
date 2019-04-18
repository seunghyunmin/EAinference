#' @title Post-selection individual inference with lasso estimator
#'
#' @description Provides confidence intervals for the set of active coefficients
#' of lasso using Metropolis-Hastings sampler.
#'
#' @param X predictor matrix.
#' @param Y response vector.
#' @param Ctype either \code{"CI"} or \code{"CS"} which represent confidence intervals
#' and confidence sets, respectively. If \code{"CI"}, confidence intervals for all active
#' coefficients are generated. If \code{"CS"}, \code{target} argument needs to be specified.
#' @param target target active variables of which one wants to generate confidence set. Needs
#' to be a subset of active set. See the example for the detail.
#' @param MHsamples optional argument. MHsamples from \code{\link{postInference.MHLS}}.
#' If MHsamples is supplied, MH sampling step is omitted. See the example for the detail.
#' @param LassoEst The result from \code{\link{lassoFit}} function with \code{type="lasso"}.
#' @param tau numeric vector. Standard deviation of proposal distribution
#'  for each beta. Adjust the value to get relevant level of acceptance rate.
#'  Default is \code{rep(1, ncol(X))}.
#' @param sig2.hat variance of error term.
#' @param alpha confidence level for confidence interval.
#' @param nChain the number of chains. For each chain, different plug-in beta will be generated
#' from its confidence region.
#' @param niterPerChain the number of iterations per chain.
#' @param method Type of robust method. Users can choose either \code{"coeff"} or \code{"mu"}.
#' @param parallel logical. If \code{parallel = TRUE}, uses parallelization.
#' Default is \code{parallel = FALSE}.
#' @param ncores integer. The number of cores to use for parallelization.
#' @param returnSamples logical. If \code{returnSamples = TRUE}, print Metropolis-Hastings samples.
#' If \code{MHsamples} is supplied, \code{returnSamples = FALSE} is forced.
#' @param ... auxiliary \code{\link{MHLS}} arguments.
#' @details
#' This function provides post-selection inference for the active coefficients selected by lasso.
#' Uses Metropolis-Hastings sampler with multiple chains to draw from the
#' distribution under a fixed active set and generates \code{(1-alpha)}
#' confidence interval for each active coefficients.
#' Set \code{returnSamples = TRUE} to check the Metropolis-Hastings samples.
#' Check the acceptance rate and adjust \code{tau} accordingly.
#' We recommend to set \code{nChain >= 10} and \code{niterPerChain >= 500}.
#'
#' @return \item{MHsamples}{a list of class MHLS.}
#' @return \item{confidenceInterval}{(1-alpha) confidence interval
#' for each active coefficient.}
#' @examples
#' set.seed(123)
#' n <- 6
#' p <- 10
#' X <- matrix(rnorm(n*p),n)
#' Y <- X %*% rep(1,p) + rnorm(n)
#' sig2 <- 1
#' lbd <- .37
#' weights <- rep(1,p)
#' LassoEst <- lassoFit(X = X, Y = Y, type = "lasso", lbd = lbd, weights = weights)
#' parallel <- (.Platform$OS.type != "windows")
#' P1 <- postInference.MHLS(LassoEst= LassoEst, X = X, Y = Y, sig2.hat = 1, alpha = .05,
#' nChain = 3, niterPerChain = 20, method = "coeff", parallel = parallel, returnSamples = TRUE)
#' P1
#' postInference.MHLS(LassoEst= LassoEst, MHsamples = P1$MHsamples,
#'                    Ctype = "CI", X = X, Y = Y, method = "coeff")
#' postInference.MHLS(LassoEst= LassoEst, MHsamples = P1$MHsamples,
#'                    Ctype = "CS", X = X, Y = Y, method = "coeff")
#' @export
postInference.MHLS <- function(LassoEst, Ctype = "CI", X, Y, sig2.hat, tau = rep(1, ncol(X)), alpha = .05,
                               MHsamples, target = which(LassoEst$B0!=0),
                               nChain = 10, method,
                               niterPerChain = 500, parallel = FALSE, ncores = 2L, returnSamples=FALSE, ...)
{
  # nChain : the number of MH chains
  # niterPerChain : the number of iteration for each chain
  # B0, S0 : The lasso estimator
  # tau : same as in MHLS function

  # LassoEst <- lassoFit(X=X, Y=Y, type = "lasso", lbd=lbd, weights=weights)
  B0 <- LassoEst$B0
  S0 <- LassoEst$S0
  lbd <- LassoEst$lbd
  weights <- LassoEst$weights

  A <- which(B0!=0)
  if (length(A)==0) {
    stop("Given lbd, active set is empty.")
  }

  if (Ctype == "CS") {
    if (any(!target %in% A)) {
      stop("target needs to be a subset of the active set A.")
    }
  }

  Y <- matrix(Y, , 1)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  nChain <- as.integer(nChain)
  niterPerChain <- as.integer(niterPerChain)

  if (missing(sig2.hat)) {
    if (length(A) >= nrow(X)) {
      stop("If size of active set matches with nrow(X), sig2.hat needs to be provided.")
    }
    sig2.hat <- summary((lm(Y~X[,A]+0)))$sigma^2
  }
  #--------------------
  # Error Handling
  #--------------------
  parallelTemp <- ErrorParallel(parallel,ncores)
  parallel <- parallelTemp$parallel
  ncores <- parallelTemp$ncores

  if (!Ctype %in% c("CI", "CS")) {
    stop("Invalid Ctype.")
  }

  if (!method %in% c("coeff", "mu")) {
    stop("Invalid method type.")
  }

  if (length(tau) != p) {
    stop("length(tau) has to be the same with col(X)")
  }

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
  if (any(c(nChain,niterPerChain) <= 0)) {
    stop("nChain & niterPerChain have to be a positive integer.")
  }

  if (length(unique(LassoEst$group))!=ncol(X)) {
    stop("lassoFit with type = 'lasso' is suppported only.")
  }

  beta.refit <- rep(0,p)
  beta.refit[A] <- coef(lm(Y~X[,A]+0))

  # if (all.equal(coef(gglasso(X, Y, pf = weights, group = 1:p, loss="ls",
  #                             intercept=F, lambda=lbd))[-1],B0) != TRUE ||
  #     all.equal(c(((t(X)/weights)%*%Y - (t(X) /weights) %*% X %*% B0) / n / lbd)
  #                , S0) != TRUE) {
  #   stop("Invalid B0 or S0, use lassoFit to get a valid lasso solution.")
  # }
  # Draw samples of pluginbeta from the 95% confidence
  #  region boundary of restricted lse.
  # If nChain ==1, we just use restricted lse.

  if (missing(MHsamples)) {
    if (method == "coeff") {
      Plugin.seq <- Pluginbeta.MHLS(X = X, Y = Y, A = A, nPlugin = nChain,
                                    sigma.hat = sqrt(sig2.hat), alpha=alpha/2)
      betaCenter <- beta.refit
      # muHat <- X[,A]%*%coef(lm(Y~X[,A]+0))
    } else {
      Plugin.seq <- PluginMu.MHLS(X = X, Y = Y, lbd = lbd, sigma.hat = sqrt(sig2.hat),
                                  # ratioSeq = seq(0,1,by=0.01), alpha = 0.05, nChain = nChain, niter = 100,
                                  alpha = alpha/2, nChain = nChain, niter = 100,
                                  method = "boundary", parallel = parallel, ncores = ncores)
      betaCenter <- rep(0,p)
      betaCenter[A] <- solve(crossprod(X[,A]))%*%t(X[,A])%*% Plugin.seq[nChain+1, ]
    }

    FF <- function(x) {
      MHLS(X = X, PE = Plugin.seq[x,], sig2 = sig2.hat, lbd = lbd,
           weights = weights, niter=niterPerChain,
           burnin = 0, B0 = B0, S0 = S0, tau = tau, PEtype = method, verbose=FALSE, ...)
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

    RefitBeta <- Refit.MHLS(X,weights,lbd,MCSAMPLE)

    if (Ctype == "CI") {
      # Using MH samples, refit the coeff.
      if (returnSamples) {
        return(list(MHsamples = TEMP, pluginValue = Plugin.seq, method = method,
                    confidenceInterval = CI.MHLS(betaRefitMH = RefitBeta, betaCenter = betaCenter,
                                                 betaRefit = beta.refit, alpha = alpha/2),
                    call = match.call()))
      } else {
        return(list(
          confidenceInterval = CI.MHLS(betaRefitMH = RefitBeta, betaRefit = beta.refit,
                                       betaCenter = betaCenter, alpha = alpha/2),
          call = match.call()))
      }
    } else {
      if (returnSamples) {
        return(list(MHsamples = TEMP, pluginValue = Plugin.seq, method = method,
                    confidenceSet = CS.MHLS(betaRefitMH = RefitBeta, betaCenter = betaCenter,
                                                 target = target, betaRefit = beta.refit, alpha = alpha/2),
                    call = match.call()))
      } else {
        return(list(
          confidenceSet = CS.MHLS(betaRefitMH = RefitBeta, betaCenter = betaCenter,
                                  target = target, betaRefit = beta.refit, alpha = alpha/2),
          call = match.call()))
      }
    }
  } else {
    if (method == "coeff") {
      betaCenter <- beta.refit
    } else {
      betaCenter <- rep(0,p)
      betaCenter[A] <- solve(crossprod(X[,A]))%*%t(X[,A])%*% Plugin.seq[nChain+1, ]
    }

    MCSAMPLE <- MHsamples[[1]]
    if (length(MHsamples) > 1) {
      for (i in 2:length(MHsamples)) {
        MCSAMPLE$beta <- rbind(MCSAMPLE$beta, MHsamples[[i]]$beta)
        MCSAMPLE$subgrad <- rbind(MCSAMPLE$subgrad, MHsamples[[i]]$subgrad)
        MCSAMPLE$acceptHistory <- MCSAMPLE$acceptHistory + MHsamples[[i]]$acceptHistory
        MCSAMPLE$niteration <- MCSAMPLE$niteration + MHsamples[[i]]$niteration
        MCSAMPLE$burnin <- MCSAMPLE$burnin + MHsamples[[i]]$burnin
      }
    }

    RefitBeta <- Refit.MHLS(X,weights,lbd,MCSAMPLE)

    if (Ctype == "CI") {
      # Using MH samples, refit the coeff.
      return(list(
        confidenceInterval = CI.MHLS(betaRefitMH = RefitBeta, betaRefit = beta.refit,
                                     betaCenter = betaCenter, alpha = alpha/2),
        call = match.call()))
    } else {
      return(list(
        confidenceSet = CS.MHLS(betaRefitMH = RefitBeta, betaCenter = betaCenter,
                                target = target, betaRefit = beta.refit, alpha = alpha/2),
        call = match.call()))

    }
  }

}


