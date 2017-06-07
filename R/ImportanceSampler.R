#' @title Computing importance weights under high-dimensional setting
#'
#' @description \code{hdIS} is used to computes importance weights using samples
#' drawn by \code{\link{DirectSampler}}. For group lasso, we provide the option
#' to use mixture distribution as a proposal distribution. See the examples
#' below for more details
#'
#' @param X Predictor matrix of size n x p, where n is the number of
#' observations and p is number of covariates.
#' @param coefftarget,sig2target,lbdtarget Parameter of target distribution.
#' (coefficient estimate, estimated variance of error, lambda)
#' @param coeffprop1,sig2prop1,lbdprop1 Parameter of proposal distribution.
#' (coefficient estimate, estimated variance of error, lambda).
#' @param coeffprop2,sig2prop2,lbdprop2 Parameter of mixture proposal
#' distribution. (coefficient estimate, estimated variance of error, lambda).
#' Only needed for group lasso under mixture distribution.
#' See Zhou and Min(2016) for more details.
#' @param propsalsample Samples drawn from proposal distribution. The samples
#' can be easily drawn by using \code{\link{DirectSampler}} with relevant choice
#' of \code{coeff, sigma2, lbd}.
#' @param group Prespecifed group structure for group lasso. p x 1 vector.
#' Use the same integer if covariates are in the same group. See example for
#' more detail.
#' @param weights Weight term for each group. Default is
#' rep(1, length(unique(group))).
#' @param TsA.method Way to construct T(eta(s),A) matrix. See Zhou and Min(2016)
#' for more detail.
#' @param log If true, importance weight is computed in log scale.
#' @param parallel logical. If true, parallelize the computation.
#' @param ncores number of CPU cores to use for parallelization.
#'
#' @details Computes importance weights which is defined as \deqn{\frac{target
#'  density}{proposal density}}, when the samples is drawn from proposal
#'  distribution with (coeffprop, sig2prop, lbdprop) while the parameter of
#'  target distribution is (coefftarget, sig2target, lbdtarget).
#'
#' @return \code{hdIS} returns importance weights of the proposed sample.
#'
#' @examples
#' set.seed(1234)
#' n <- 10
#' p <- 30
#' Niter <-  10
#' Group <- rep(1:(p/10), each = 10)
#' Weights <- rep(1, p/10)
#' x <- matrix(rnorm(n*p), n)
#'
#' # Target distribution parameter
#' coefftarget <- rep(0,p)
#' sig2target <- .5
#' lbdtarget <- .37
#'
#' # Proposal distribution parameter
#' coeffprop1 <- rep(0,p) ;coeffprop2 <- rep(1,p)
#' sig2prop1 <- .5; sig2prop2 <- 1
#' lbdprop1 <- .37; lbdprop2 <- .5
#'
#' #
#' # Using non-mixture distribution
#' # ------------------------------
#' # Target distribution parameters (coeff, sig2, lbd) = (rep(0,p), .5, .37)
#' # Proposal distribution parameters (coeff, sig2, lbd) = (rep(1,p), 1, .5)
#' #
#' DS <- DirectSampler(X = x, coeff_1 = rep(1, p), sig2_1 = 1, lbd_1 = .5,
#'  weights = Weights, group = Group, niter = Niter, parallel = FALSE)
#'
#' hdIS(X = x,coefftarget = rep(0,p), sig2target = .5, lbdtarget = .37,
#'      coeffprop1 = rep(1,p),sig2prop1 = 1,lbdprop1 = .5,proposalsample = DS, group = Group,
#'      weights = Weights, log = TRUE)
#'
#' #
#' # Using mixture distribution
#' # ------------------------------
#' # Target distribution parameters (coeff, sig2, lbd) = (rep(0,p), .5, .37)
#' # Proposal distribution parameters
#' #  (coeff, sig2, lbd) = (rep(0,p), .5, .37) & (rep(1,p), 1, .5)
#' #
#' #
#' coefftarget <- rep(0,p)
#' sig2target <- .5
#' lbdtarget <- .37
#' coeffprop1 <- rep(0,p)
#' coeffprop2 <- rep(1,p)
#' sig2prop1 <- .5
#' sig2prop2 <- 1
#' lbdprop1 <- .37
#' lbdprop2 <- .5
#'
#' DSMixture <- DirectSampler(X = x, coeff_1 = coeffprop1, sig2_1 = sig2prop1, lbd_1 = lbdprop1,
#'  coeff_2 = coeffprop2, sig2_2 = sig2prop2, lbd_2 = lbdprop2, weights = Weights,
#'  group = Group, niter = Niter, parallel = TRUE)
#' hdIS(X = x,coefftarget = coefftarget, sig2target = sig2target, lbdtarget = lbdtarget,
#'      coeffprop1 = coeffprop1, sig2prop1 = sig2prop1, lbdprop1 = lbdprop1,
#'      coeffprop2 = coeffprop2, sig2prop2 = sig2prop2, lbdprop2 = lbdprop2,
#'      proposalsample = DSMixture, group = Group,
#'      weights = Weights, log = TRUE)
#' @export
hdIS=function(X, coefftarget, sig2target, lbdtarget, coeffprop1, sig2prop1,
  lbdprop1, coeffprop2, sig2prop2, lbdprop2, proposalsample, group,
  weights = rep(1, length(unique(group))), TsA.method = "default", log = TRUE,
  parallel = FALSE, ncores)
{
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (parallel && !missing(ncores) && ncores == 1) {
    ncores <- 2
    warning("If parallel=TRUE, ncores needs to be greater than 1. Automatically
            Set ncores to the maximum number.")
  }

  if (any(c(length(coefftarget), length(coeffprop1)) != p) ||
      (!missing(coeffprop2) && (length(coeffprop2) != p)) ) {
    stop("coefftarget/coeffprop must have a same length with the number of X columns")
  }

  if (length(group) != p) {
    stop("group must have a same length with the number of X columns")
  }
  if (length(weights) != length(unique(group))) {
    stop("weights has to have a same length as the number of groups")
  }

  if (all(group==1:p)) {
    #
    # Lasso
    #
    # precalculation
    C <- t(X) %*% X / n
    egC <- eigen(C)
    V <- egC$vectors
    R <- 1:n
    N <- (n+1):p
    InvVarR     <- 1 / (egC$values[R] * sig2target / n) #inverse of (sig2target*Lambda_i/n)
    InvVarRprop <- 1 / (egC$values[R] * sig2prop1 / n) #inverse of (sig2prop1*Lambda_i/n)
    VR <-matrix(V[, R], p, n)
    VRC <- VRCprop <- t(VR)%*%C
    W <- diag(weights)
    LBD <- LBDprop <- diag(egC$values[R])
    VRW <- VRWprop <- t(VR)%*%W
    VRCB     <- t(VR) %*% C %*% coefftarget
    VRCBprop <- t(VR) %*% C %*% coeffprop1

    #sample from proposal distribution
    B <- proposalsample$beta
    S <- proposalsample$subgrad
    niter <- nrow(B)

    logISweights <- numeric(niter)
    FF <- function(x) {
      (n - sum(B[x,] != 0)) * log(lbdtarget / lbdprop1) - 0.5 * sum((VRC %*% B[x, ] + lbdtarget * VRW %*% S[x, ] -
      VRCB)^2 * InvVarR) + 0.5 * sum((VRCprop %*% B[x, ] + lbdprop1 * VRWprop %*% S[x, ] - VRCBprop)^2 * InvVarRprop)
    }
    if (!parallel) {
      for (t in 1:niter) {
        logISweights[t] <- FF(t)
      }
    } else {
      if (missing(ncores)) {
        ncores <- parallel::detectCores()
      } else if (ncores > parallel::detectCores()){
        ncores <- parallel::detectCores()
        warnings("ncores is larger than the maximum number of available processes.
                 Set it to the maximum possible value.")
      }
      logISweights <- parallel::mclapply(1:niter, FF, mc.cores = ncores)
      logISweights <- do.call(c,Weight)
     }
    logISweights <- logISweights - n / 2 * (log(sig2target / sig2prop1))
    return(ifelse(log, logISweights, exp(logISweights)))
  } else {
    #
    # Group Lasso
    #
    if (!TsA.method %in% c("default", "qr")) {
      stop("TsA.method should be either \"default\" or \"qr\"")
    }
    # if (length(sig2prop) != length(lbdprop) || length(sig2prop) != nrow(coeffprop) || length(lbdprop) != nrow(coeffprop)) {
    #   stop("provide all the parameters for the proposal distribution(s)")
    # }
    #
    if (!sum(c(missing(sig2prop2), missing(lbdprop2), missing(coeffprop2)))
        %in% c(0,3)) {
      stop("provide all the parameters for the proposal distribution(s)")
    }

    if (missing(sig2prop2)) {
      Mixture <- FALSE
      lbdprop2 <- lbdprop1
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

    if (!all(lbdtarget == c(lbdprop1, lbdprop2)) && TsA.method != "null") {
      Q <- MASS::Null(t(X)/W)#round(Null(t(X)/W),5)
    }
    t.XWinv <- t(X)/W
    Weight <- c()

    FF <- function(x) {
      Beta <- B[x,]
      Subgrad <- S[x,]
      if (n < p) {
        H.tilde.target <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - coefftarget) + lbdtarget * W * Subgrad) #H.tilde
        H.tilde.prop1 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - coeffprop1) + lbdprop1 * W * Subgrad) #H.tilde proposed1
        if (Mixture) {
          H.tilde.prop2 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - coeffprop2) + lbdprop2 * W * Subgrad) #H.tilde proposed2
        }

        r <- group.norm2(Beta, group)
        A <- unique(group[Beta != 0])

        if (!all(lbdtarget == c(lbdprop1, lbdprop2))) {

          if (TsA.method == "null") {
            TSA <- TsA.select(t.XWinv, Subgrad, group, A, n, p)
          } else {
            TSA <- TsA.select(Q, Subgrad, group, A, n, p)
          }

          log.f1 <- sum(dnorm(H.tilde.prop1, 0, sqrt(sig2prop1/n), log = T)) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdprop1, weights, TSA) )
          if (Mixture) {
            log.f2 <- sum(dnorm(H.tilde.prop2, 0, sqrt(sig2prop2/n), log = T)) +
              (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdprop2, weights, TSA) )
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
          log.f1 <- sum(dnorm(H.tilde.prop1, 0, sqrt(sig2prop1/n), log = TRUE))
          if (Mixture) {log.f2 <- sum(dnorm(H.tilde.prop2, 0, sqrt(sig2prop2/n), log = TRUE))}
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
        H.prop1 <-  (Psi %*% (Beta - coeffprop1) + lbdprop1 * W * Subgrad) #H.tilde proposed1
        if (Mixture) H.prop2 <-  (Psi %*% (Beta - coeffprop2) + lbdprop2 * W * Subgrad) #H.tilde proposed2

        r <- group.norm2(Beta, group)
        A <- unique(group[Beta != 0])

        if (!all(lbdtarget == c(lbdprop1, lbdprop2))) {

          if (TsA.method == "null") {
            TSA <- TsA.select(t.XWinv, Subgrad, group, A, n, p)
          } else {
            TSA <- TsA.select(Q, Subgrad, group, A, n, p)
          }

          log.f1 <- dmvnorm(drop(H.prop1), , sig2prop1/n * Psi, log = T) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdprop1, weights, TSA) )
          if (Mixture) {
            log.f2 <- dmvnorm(drop(H.prop2), , sig2prop2/n * Psi, log = T) +
              (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdprop2, weights, TSA) )
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
          log.f1 <- sum(dmvnorm(drop(H.prop1), , sig2prop1/n * Psi, log = TRUE))
          if (Mixture) {
            log.f2 <- sum(dmvnorm(drop(H.prop2), , sig2prop2/n * Psi, log = TRUE))
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
      return(Weight)
    }

    if (!parallel) {
      for (t in 1:niter) {
        Weight[t] <- FF(t)
      }
    } else {
      if (missing(ncores)) {
        ncores <- parallel::detectCores()
      } else if (ncores > parallel::detectCores()){
        ncores <- parallel::detectCores()
        warnings("ncores is larger than the maximum number of available processes.
                 Set it to the maximum possible value.")
      }
      Weight <- parallel::mclapply(1:niter, FF, mc.cores = ncores)
      Weight <- do.call(c,Weight)
    }
    return(Weight)
  }
}
