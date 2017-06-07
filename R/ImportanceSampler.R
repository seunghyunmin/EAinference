#' With samples drawn from proposal distribution, this function computes the importance weight.
#' Simulating under proposal distribution can be easily done by using EALasso::DirectSampler.
#'
#' @param X Predictor matrix of size n x p, where n is the number of observations and p is number of covariates.
#' @param coefftarget Parameter of target distribution. n x 1 vector of coefficient estimate.
#' @param sig2target Parameter of target distribution. Estimated variance of error term.
#' @param lbdtarget Parameter of target distribution. Lambda.
#' @param coeffprop Coefficient estimate of proposal distribution.
#' @param sig2prop Estimated variance of error term of proposal distribution.
#' @param lbdprop Lambda of proposal distribution. Lambda.
#' @param propsalsample Samples drawn from proposal distribution. The samples can be easily drawn by using DirectSampler.
#' @param group Prespecifed group structure for group lasso. p x 1 vector.
#' Use the same integer if covariates are in the same group. See example for more detail.
#' @param weights Weight term for each group. Default is rep(1, length(unique(group))).
#' @param TsA.method Way to construct T(eta(s),A) matrix. See Zhou and Min(2016) for more detail.
#' @param log If true, importance weight is computed in log scale.
#' @param parallel If true, parallelize the computation.
#' @param ncores number of CPU cores to use for parallelization.
#'
#' @details If futype="normal", it generate
#' @return \describe{
#'   \item{beta}{coefficient matrix of size N x p.}
#'   \item{subgrad}{subgradient matrix of size N x p.}
#'  }
#' @examples
#' n <- 10
#' p <- 30
#' lbd <-  .5
#' niter <-  10
#' group <- rep(1:(p/10), each=10)
#' Weights <- rep(1,p/10)
#' x <- matrix(rnorm(n*p), n)
#' DirectSampler(X = x, coefftarget = rep(0,p), sig2target=1, lbd=.5, weights=Weights, group=group,N=niter, parallel=FALSE)
#' DirectSampler(X = x, coefftarget = rep(0,p), sig2target=1, lbd=.5, weights=Weights, group=group,N=niter, parallel=TRUE)
#' @export
hdIS=function(X, coefftarget, sig2target, lbdtarget, coeffprop, sig2prop, lbdprop, proposalsample,
   group, weights=rep(1,length(unique(group))), TsA.method = "default", log=TRUE, parallel=FALSE, ncores)
{
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (length(coefftarget) != p || length(coeffprop)%%p !=0) {
    stop("coefftarget/coeffprop must have a same length with the number of X columns")
  }
  if (length(group) != p) {
    stop("group must have a same length with the number of X columns")
  }
  if (length(weights) != length(unique(group))) {
    stop("weights has to have a same length as the number of groups")
  }

  if (all(group==1:p)) {
    # Lasso
    # precalculation
    C <- t(X) %*% X / n
    egC <- eigen(C)
    V <- egC$vectors
    R <- 1:n
    N <- (n+1):p
    InvVarR     <- 1 / (egC$values[R] * sig2target / n) #inverse of (sig2target*Lambda_i/n)
    InvVarRprop <- 1 / (egC$values[R] * sig2prop / n) #inverse of (sig2prop*Lambda_i/n)
    VR <-matrix(V[, R], p, n)
    VRC <- VRCprop <- t(VR)%*%C
    W <- diag(weights)
    LBD <- LBDprop <- diag(egC$values[R])
    VRW <- VRWprop <- t(VR)%*%W
    VRCB     <- t(VR) %*% C %*% coefftarget
    VRCBprop <- t(VR) %*% C %*% coeffprop

    #sample from proposal distribution
    B <- proposalsample$beta
    S <- proposalsample$subgrad
    niter <- nrow(B)

    logISweights <- numeric(niter)
    # for(t in 1:niter)
    # {
    #   logISweights[t] <- (n-sum(B[t,]!=0))*log(lbdtarget/lbdprop)-0.5*sum((VRC%*%B[t,]+lbdtarget*VRW%*%S[t,]-VRCB)^2*InvVarR)+
    #     0.5*sum((VRCprop%*%B[t,]+lbdprop*VRWprop%*%S[t,]-VRCBprop)^2*InvVarRprop)
    # }
    # logISweights <- logISweights-n/2*(log(sig2target/sig2prop))
    FF <- function(x) {
      (n - sum(B[x,] != 0)) * log(lbdtarget / lbdprop) - 0.5 * sum((VRC %*% B[x, ]+lbdtarget * VRW %*% S[x, ] -
      VRCB)^2 * InvVarR) + 0.5 * sum((VRCprop %*% B[x, ] + lbdprop * VRWprop %*% S[x, ] - VRCBprop)^2 * InvVarRprop)
    }
    if (!parallel) {
      for (t in 1:niter) {
        logISweights[t] <- FF(t)
      }
    } else {
      if (missing(ncores)) {ncores <- detectCores()}
      options(mc.cores <- ncores)
      logISweights <- mclapply(1:niter, FF)
      logISweights <- do.call(c, Weight)
    }
    logISweights <- logISweights - n / 2 * (log(sig2target / sig2prop))
    return(if(log) logISweights else exp(logISweights))
  } else {
    # Group Lasso
    coeffprop <- matrix(coeffprop, , p)
    if (!TsA.method %in% c("default", "qr")) {
      stop("TsA.method should be either \"default\" or \"qr\"")
    }
    if (length(sig2prop) != length(lbdprop) || length(sig2prop) != nrow(coeffprop) || length(lbdprop) != nrow(coeffprop)) {
      stop("provide all the parameters for the proposal distribution(s)")
    }

    if (length(sig2prop) == 1) {
      Mixture <- FALSE
      lbdprop[2] <- lbdprop[1]
    } else {
      Mixture <- TRUE
    }
    #TsA.select <- switch(TsA.method, null = TsA.null, qr = TsA.qr, default = TsA)
    TsA.select <- switch(TsA.method, qr = TsA.qr, default = TsA)
    Psi <- t(X)%*%X / n
    S <- proposalsample$subgrad
    B <- proposalsample$beta
    niter <- nrow(B)

    #precaculation
    if (n < p) {
      ginv.tX <- solve(tcrossprod(X)) %*% X
    }

    W <- rep(weights, table(group))

    if (!all(lbdtarget == c(lbdprop[1], lbdprop[2])) && TsA.method != "null") {
      Q <- Null(t(X)/W)#round(Null(t(X)/W),5)
    }
    t.XWinv <- t(X)/W
    Weight <- c()

    FF <- function(x) {
      Beta <- B[x,]
      Subgrad <- S[x,]
      if (n < p) {
        H.tilde.target <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - coefftarget) + lbdtarget * W * Subgrad) #H.tilde
        H.tilde.prop1 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - coeffprop[1,]) + lbdprop[1] * W * Subgrad) #H.tilde proposed1
        if (Mixture) {
          H.tilde.prop2 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - coeffprop[2,]) + lbdprop[2] * W * Subgrad) #H.tilde proposed2
        }

        r <- group.norm2(Beta, group)
        A <- unique(group[Beta != 0])

        if (!all(lbdtarget == c(lbdprop[1], lbdprop[2]))) {

          if (TsA.method == "null") {
            TSA <- TsA.select(t.XWinv, Subgrad, group, A, n, p)
          } else {
            TSA <- TsA.select(Q, Subgrad, group, A, n, p)
          }

          log.f1 <- sum(dnorm(H.tilde.prop1, 0, sqrt(sig2prop[1]/n), log = T)) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdprop[1], weights, TSA) )
          if (Mixture) {
            log.f2 <- sum(dnorm(H.tilde.prop2, 0, sqrt(sig2prop[2]/n), log = T)) +
              (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdprop[2], weights, TSA) )
          }
          log.f0 <- sum(dnorm(H.tilde.target, 0, sqrt(sig2target/n), log = T)) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdtarget, weights, TSA) )

          if (!Mixture) {
            if (log) {Weight <- log.f0 - log.f1} else {
              Weight <- exp(log.f0 - log.f1)
            }
          } else {
            if (log) {
              Weight <- - log(exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
            } else {
              Weight <- 1 / (exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
            }
          }
        } else {
          log.f1 <- sum(dnorm(H.tilde.prop1, 0, sqrt(sig2prop[1]/n), log = TRUE))
          if (Mixture) {log.f2 <- sum(dnorm(H.tilde.prop2, 0, sqrt(sig2prop[2]/n), log = TRUE))}
          log.f0 <- sum(dnorm(H.tilde.target, 0, sqrt(sig2target/n), log = TRUE))

          if (!Mixture) {
            if (log) {Weight <- log.f0 - log.f1} else {
              Weight <- exp(log.f0 - log.f1)
            }
          } else {
            if (log) {
              Weight <- - log(exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
            } else {
              Weight <- 1/ (exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
            }
          }
        }
      }
      else {
        H.target <- (Psi %*% (Beta - coefftarget) + lbdtarget * W * Subgrad) #H.tilde
        H.prop1 <-  (Psi %*% (Beta - coeffprop[1,]) + lbdprop[1] * W * Subgrad) #H.tilde proposed1
        if (Mixture) H.prop2 <-  (Psi %*% (Beta - coeffprop[2,]) + lbdprop[2] * W * Subgrad) #H.tilde proposed2

        r <- group.norm2(Beta, group)
        A <- unique(group[Beta != 0])

        if (!all(lbdtarget == c(lbdprop[1], lbdprop[2]))) {

          if (TsA.method == "null") {
            TSA <- TsA.select(t.XWinv, Subgrad, group, A, n, p)
          } else {
            TSA <- TsA.select(Q, Subgrad, group, A, n, p)
          }

          log.f1 <- dmvnorm(drop(H.prop1), , sig2prop[1]/n * Psi, log = T) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdprop[1], weights, TSA) )
          if (Mixture) {
            log.f2 <- dmvnorm(drop(H.prop2), , sig2prop[2]/n * Psi, log = T) +
              (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdprop[2], weights, TSA) )
          }
          log.f0 <- dmvnorm(drop(H.target), , sig2target/n * Psi, log = T) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdtarget, weights, TSA) )

          if (!Mixture) {
            if (log) {Weight <- log.f0 - log.f1} else {
              Weight <- exp(log.f0 - log.f1)
            }
          } else {
            if (log) {
              Weight <- - log(exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
            } else {
              Weight <- 1/ (exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
            }
          }
        } else {
          log.f1 <- sum(dmvnorm(drop(H.prop1), , sig2prop[1]/n * Psi, log = TRUE))
          if (Mixture) {
            log.f2 <- sum(dmvnorm(drop(H.prop2), , sig2prop[2]/n * Psi, log = TRUE))
          }
          log.f0 <- sum(dmvnorm(drop(H.target), , sig2target/n * Psi, log = TRUE))

          if (!Mixture) {
            if (log) {Weight <- log.f0 - log.f1} else {
              Weight <- exp(log.f0 - log.f1)
            }
          } else {
            if (log) {
              Weight <- - log(exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
            } else {
              Weight <- 1/ (exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
            }
          }
        }
      }
    }

    if (!parallel) {
      for (t in 1:niter) {
        Weight[t] <- FF(t)
      }
    } else {
      if (missing(ncores)) {ncores <- detectCores()}
      options(mc.cores=ncores)
      Weight <- mclapply(1:niter,FF)
      Weight <- do.call(c,Weight)
    }
    return(Weight)
  }
}
