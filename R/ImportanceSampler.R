#' @title Computing importance weights under high-dimensional setting
#'
#' @description \code{hdIS} is used to computes importance weights using samples
#' drawn by \code{\link{DirectSampler}}. For group lasso, we provide the option
#' to use mixture distribution as a proposal distribution. See the examples
#' below for more details
#'
#' @param DirectSample Bootstrap samples of class \code{DS} from \code{DirectSampler}
#' @param pluginTarget,sig2Target,lbdTarget Parameter of target distribution.
#' (coefficient estimate, estimated variance of error, lambda)
#' @param TsA.method Way to construct T(eta(s),A) matrix. See Zhou and Min(2016)
#' for more detail.
#' @param log If true, importance weight is computed in log scale.
#' @param parallel logical. If true, parallelize the computation.
#' @param ncores number of CPU cores to use for parallelization.
#'
#' @details Computes importance weights which is defined as \deqn{\frac{target
#'  density}{proposal density}}, when the samples is drawn from proposal
#'  distribution with (coeffprop, sig2prop, lbdprop) while the parameter of
#'  target distribution is (pluginTarget, sig2Target, lbdTarget).
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
#' pluginTarget <- rep(0,p)
#' sig2Target <- .5
#' lbdTarget <- .37
#'
#' # Proposal distribution parameter
#' pluginProp1 <- rep(0,p) ;pluginProp2 <- rep(1,p)
#' sig2Prop1 <- .5; sig2Prop2 <- 1
#' lbdProp1 <- .37; lbdProp2 <- .5
#'
#' #
#' # Using non-mixture distribution
#' # ------------------------------
#' # Target distribution parameters (coeff, sig2, lbd) = (rep(0,p), .5, .37)
#' # Proposal distribution parameters (coeff, sig2, lbd) = (rep(1,p), 1, .5)
#' #
#' DS <- DirectSampler(X = x, pointEstimate_1 = rep(1, p), sig2_1 = 1, lbd_1 = .5,
#'  weights = Weights, group = Group, niter = Niter, type = "coeff", parallel = FALSE)
#'
#' hdIS(DS, pluginTarget = rep(0,p), sig2Target = .5, lbdTarget = .37,
#'  log = TRUE)
#'
#' #
#' # Using mixture distribution
#' # ------------------------------
#' # Target distribution parameters (coeff, sig2, lbd) = (rep(0,p), .5, .37)
#' # Proposal distribution parameters
#' #  (coeff, sig2, lbd) = (rep(0,p), .5, .37) & (rep(1,p), 1, .5)
#' #
#' #
#' pluginProp1 <- rep(0,p)
#' pluginProp2 <- rep(1,p)
#' sig2Prop1 <- .5
#' sig2Prop2 <- 1
#' lbdProp1 <- .37
#' lbdProp2 <- .5
#'
#' DSMixture <- DirectSampler(X = x, pointEstimate_1 = pluginProp1,
#'  sig2_1 = sig2Prop1, lbd_1 = lbdProp1, pointEstimate_2 = pluginProp2,
#'  sig2_2 = sig2Prop2, lbd_2 = lbdProp2, weights = Weights, group = Group,
#'  niter = Niter, type = "coeff", parallel = TRUE)
#' hdIS(DSMixturem, pluginTarget = rep(0,p), sig2Target = .5, lbdTarget = .37,
#'  log = TRUE)
#' @export
hdIS=function(DirectSample, pluginTarget, sig2Target, lbdTarget,
            TsA.method = "default", log = TRUE, parallel = FALSE, ncores = 2L)
{
  if (class(DirectSample) != "DS") {
    stop("Use EAlasso::DirectSampler to generate DirectSample.")
  }

  X <- DirectSample$X
  n <- nrow(X)
  p <- ncol(X)

  if (DirectSample$mixture) {
    pluginProp1 <- DirectSample$PointEstimate[1,]
    pluginProp2 <- DirectSample$PointEstimate[2,]
    sig2Prop1 <- DirectSample$sig2[1]
    sig2Prop2 <- DirectSample$sig2[2]
    lbdProp1 <- DirectSample$lbd[1]
    lbdProp2 <- DirectSample$lbd[2]
  } else {
    pluginProp1 <- DirectSample$PointEstimate
    sig2Prop1 <- DirectSample$sig2
    lbdProp1 <- DirectSample$lbd
  }

  type <- DirectSample$type
  method <- DirectSample$method
  group <- DirectSample$group
  weights <- DirectSample$weights

  B <- DirectSample$beta
  S <- DirectSample$subgrad
  niter <- nrow(B)

  if (n >= p) {
    stop("High dimensional setting is required, i.e. nrow(X) < ncol(X) required.")
  }

  if (parallel && ncores == 1) {
    ncores <- 2
    warning("If parallel=TRUE, ncores needs to be greater than 1. Automatically
            Set ncores to the maximum number.")
  }
  if (parallel && ncores > parallel::detectCores()){
    ncores <- parallel::detectCores()
    warning("ncores is larger than the maximum number of available processes.
             Set it to the maximum possible value.")
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
    InvVarR     <- 1 / (egC$values[R] * sig2Target / n) #inverse of (sig2Target*Lambda_i/n)
    InvVarRprop <- 1 / (egC$values[R] * sig2Prop1 / n) #inverse of (sig2Prop1*Lambda_i/n)
    VR <-matrix(V[, R], p, n)
    VRC <- VRCprop <- t(VR)%*%C
    W <- diag(weights)
    LBD <- LBDprop <- diag(egC$values[R])
    VRW <- VRWprop <- t(VR)%*%W
    if (type == "coeff") {
      VRCB     <- t(VR) %*% C %*% pluginTarget
      VRCBprop <- t(VR) %*% C %*% pluginProp1
    } else {
      VRCB     <- t(VR) %*% t(X) %*% pluginTarget / n
      VRCBprop <- t(VR) %*% t(X) %*% pluginProp1 / n
    }

    #sample from proposal distribution

    logISweights <- numeric(niter)
    FF <- function(x) {
      (n - sum(B[x,] != 0)) * log(lbdTarget / lbdProp1) - 0.5 * sum((VRC %*% B[x, ] + lbdTarget * VRW %*% S[x, ] -
      VRCB)^2 * InvVarR) + 0.5 * sum((VRCprop %*% B[x, ] + lbdProp1 * VRWprop %*% S[x, ] - VRCBprop)^2 * InvVarRprop)
    }
    if (!parallel) {
      for (t in 1:niter) {
        logISweights[t] <- FF(t)
      }
    } else {
      logISweights <- parallel::mclapply(1:niter, FF, mc.cores = ncores)
      logISweights <- do.call(c,Weight)
     }
    logISweights <- logISweights - n / 2 * (log(sig2Target / sig2Prop1))
    return(ifelse(log, logISweights, exp(logISweights)))
  } else {
    #
    # Group Lasso
    #
    if (!TsA.method %in% c("default", "qr")) {
      stop("TsA.method should be either \"default\" or \"qr\"")
    }

    if (!DirectSample$mixture) {
      Mixture <- FALSE
      lbdProp2 <- lbdProp1
    } else {
      Mixture <- TRUE
    }
    #TsA.select <- switch(TsA.method, null = TsA.null, qr = TsA.qr, default = TsA)
    TsA.select <- switch(TsA.method, qr = TsA.qr, default = TsA)
    Psi <- t(X)%*%X / n

    #precaculation
    if (n < p) {
      ginv.tX <- solve(tcrossprod(X)) %*% X
    }

    W <- rep(weights, table(group))

    if (!all(lbdTarget == c(lbdProp1, lbdProp2)) && TsA.method != "null") {
      Q <- MASS::Null(t(X)/W)#round(Null(t(X)/W),5)
    }
    t.XWinv <- t(X)/W
    Weight <- numeric(niter)

    FF <- function(x) {
      Beta <- B[x,]
      Subgrad <- S[x,]
      if (n < p) {
        if (type == "coeff") {
          H.tilde.target <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - pluginTarget) + lbdTarget * W * Subgrad) #H.tilde
          H.tilde.prop1 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - pluginProp1) + lbdProp1 * W * Subgrad) #H.tilde proposed1
          if (Mixture) {
            H.tilde.prop2 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - pluginProp2) + lbdProp2 * W * Subgrad) #H.tilde proposed2
          }
        } else {
          H.tilde.target <- sqrt(n) * ginv.tX %*% (Psi %*% Beta + lbdTarget * W * Subgrad) - pluginTarget / sqrt(n) #H.tilde
          H.tilde.prop1 <- sqrt(n) * ginv.tX %*% (Psi %*% Beta + lbdProp1 * W * Subgrad) - pluginProp1 / sqrt(n)  #H.tilde proposed1
          if (Mixture) {
            H.tilde.prop2 <- sqrt(n) * ginv.tX %*% (Psi %*% Beta + lbdProp2 * W * Subgrad) - pluginProp2 / sqrt(n)   #H.tilde proposed2
          }
        }

        r <- group.norm2(Beta, group)
        A <- unique(group[Beta != 0])

        if (!all(lbdTarget == c(lbdProp1, lbdProp2))) {
#           if (TsA.method == "null") {
#             TSA <- TsA.select(t.XWinv, Subgrad, group, A, n, p)
#           } else {
#             TSA <- TsA.select(Q, Subgrad, group, A, n, p)
#           }
          TSA <- TsA.select(Q, Subgrad, group, A, n, p)

          log.f1 <- sum(dnorm(H.tilde.prop1, 0, sqrt(sig2Prop1/n), log = T)) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdProp1, weights, TSA) )
          if (Mixture) {
            log.f2 <- sum(dnorm(H.tilde.prop2, 0, sqrt(sig2Prop2/n), log = T)) +
              (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdProp2, weights, TSA) )
          }
          log.f0 <- sum(dnorm(H.tilde.target, 0, sqrt(sig2Target/n), log = T)) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdTarget, weights, TSA) )

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
          log.f1 <- sum(dnorm(H.tilde.prop1, 0, sqrt(sig2Prop1/n), log = TRUE))
          if (Mixture) {log.f2 <- sum(dnorm(H.tilde.prop2, 0, sqrt(sig2Prop2/n), log = TRUE))}
          log.f0 <- sum(dnorm(H.tilde.target, 0, sqrt(sig2Target/n), log = TRUE))

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
        if (type == "coeff") {
          H.target <- Psi %*% (Beta - pluginTarget) + lbdTarget * W * Subgrad #H.tilde
          H.prop1 <-  Psi %*% (Beta - pluginProp1) + lbdProp1 * W * Subgrad #H.tilde proposed1
          if (Mixture) H.prop2 <-  Psi %*% (Beta - pluginProp2) + lbdProp2 * W * Subgrad #H.tilde proposed2
        } else {
          H.target <- Psi %*% Beta + lbdTarget * W * Subgrad - t(X) %*% pluginTarget / n #H.tilde
          H.prop1 <-  Psi %*% Beta + lbdProp1 * W * Subgrad - t(X) %*% pluginProp1 / n  #H.tilde proposed1
          if (Mixture) H.prop2 <-  Psi %*% Beta + lbdProp2 * W * Subgrad - t(X) %*% pluginProp2 / n  #H.tilde proposed2
        }

        r <- group.norm2(Beta, group)
        A <- unique(group[Beta != 0])

        if (!all(lbdTarget == c(lbdProp1, lbdProp2))) {

          if (TsA.method == "null") {
            TSA <- TsA.select(t.XWinv, Subgrad, group, A, n, p)
          } else {
            TSA <- TsA.select(Q, Subgrad, group, A, n, p)
          }

          log.f1 <- dmvnorm(drop(H.prop1), , sig2Prop1/n * Psi, log = T) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdProp1, weights, TSA) )
          if (Mixture) {
            log.f2 <- dmvnorm(drop(H.prop2), , sig2Prop2/n * Psi, log = T) +
              (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdProp2, weights, TSA) )
          }
          log.f0 <- dmvnorm(drop(H.target), , sig2Target/n * Psi, log = T) +
            (log.Jacobi.partial(X, Subgrad, r, Psi, group, A, lbdTarget, weights, TSA) )

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
          log.f1 <- sum(dmvnorm(drop(H.prop1), , sig2Prop1/n * Psi, log = TRUE))
          if (Mixture) {
            log.f2 <- sum(dmvnorm(drop(H.prop2), , sig2Prop2/n * Psi, log = TRUE))
          }
          log.f0 <- sum(dmvnorm(drop(H.target), , sig2Target/n * Psi, log = TRUE))

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
      Weight <- parallel::mclapply(1:niter, FF, mc.cores = ncores)
      Weight <- do.call(c,Weight)
    }
    return(Weight)
  }
}
# hdIS=function(X, pluginTarget, sig2Target, lbdTarget, pluginProp1, sig2Prop1,
#   lbdProp1, pluginProp2, sig2Prop2, lbdProp2, proposalsample, group,
#   weights = rep(1, length(unique(group))), type = "coeff", TsA.method = "default", log = TRUE,
#   parallel = FALSE, ncores)
