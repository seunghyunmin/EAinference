#' @title Parameteric Bootstrapping lasso / group lasso estimator
#'
#' @description Draw bootstrap samples in parametric way and
#' derive (group) lasso estimator along with subgradient.
#'
#' @param X Predictor matrix.
#' @param pointEstimate_1,sig2_1,lbd_1 Parameter of target distribution.
#' (Estimate of true coefficient or E(y) depends on \code{PEtype}, estimated variance of error, lambda).
#' @param pointEstimate_2,sig2_2,lbd_2 Additional Parameter of target distribution only
#' if mixture distribution is used.
#' @param weights Weight vector in length of the number of groups. Default is
#' \code{rep(1, max(group))}.
#' @param group p x 1 vector of consecutive integers describing group structure.
#' The number of groups should be same as max(group). Default is \code{group = 1:p}
#' , where \code{p} is number of covariates.
#' @param niter The number of iterations.
#' @param type type of penalty, either to be "lasso", "grlasso", "slasso" or "sgrlasso".
#' @param PEtype either to be "coeff" or "mu". Decide what kind of \code{pointEstimate} to use.
#' @param parallel Logical. If \code{TRUE}, use parallelization.
#' @param ncores Integer. The number of cores to use for the parallelization.
#' @param verbose Whether to show the process. Default is FALSE. Only works when
#'  parallel=FALSE.
#'
#' @details This function provides bootstrap samples for (group) lasso estimator
#' and its subgradient. The sampling distribution is chracterized by \code{(pointEstimate, sig2, lbd)}.
#' First, we generate \code{y_new} by \code{X * pointEstimate + error_new}, while \code{error_new}
#' is generated from N(0, sig2). See Zhou(2014) and Zhou and Min(2016) for more details.
#'
#' Four distict penalties can be used; "lasso" for lasso, "grlasso" for group lasso, "slasso" for scaled lasso
#' and "sgrlasso" for scaled group lasso.
#'
#' If non-mixture distribution is used, the distribution with parameters \code{(pointEstimate_1, sig2_1, lbd_1)}
#' will be used.
#' If one uses mixture distribution by providing \code{(pointEstimate_2, sig2_2, lbd_2)},
#' with 1/2 chance, samples will be drawn from the distribution with
#' (pointEstimate_1, sig2_1, lbd_1) and with another 1/2 chance, they will be drawn from
#' the distribution with (pointEstimate_2, sig2_2, lbd_2).
#'
#' @return \item{beta}{(group) lasso estimator.}
#' @return \item{subgrad}{subgradient.}
#' @return \item{hsigma}{standard deviation estimator, for type="slasso" or type="sgrlasso" only.}
#' @return \item{X, pointEstimate, sig2, weights, group, type, PEtype, mixture}{Model parameters}
#' @examples
#' set.seed(1234)
#' n <- 10
#' p <- 30
#' Niter <-  10
#' Group <- rep(1:(p/10), each = 10)
#' Weights <- rep(1, p/10)
#' x <- matrix(rnorm(n*p), n)
#' #
#' # Using non-mixture distribution
#' #
#' PBsampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
#'  weights = Weights, group = Group, type = "grlasso", niter = Niter, parallel = FALSE)
#' PBsampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
#'  weights = Weights, group = Group, type = "grlasso", niter = Niter, parallel = TRUE)
#' #
#' # Using mixture distribution
#' #
#' PBsampler(X = x, pointEstimate_1 = rep(0, p), sig2_1 = 1, lbd_1 = .5,
#'  pointEstimate_2 = rep(1, p), sig2_2 = 2, lbd_2 = .3, weights = Weights,
#'  group = Group, type = "grlasso", niter = Niter, parallel = TRUE)
#' @export
PBsampler <- function(X, pointEstimate_1, sig2_1, lbd_1, pointEstimate_2,
  sig2_2, lbd_2, weights = rep(1, max(group)), group = 1:ncol(X), niter = 2000,
  type, PEtype = "coeff", parallel = FALSE, ncores = 2L,
  verbose = FALSE)
{
  n <- nrow(X)
  p <- ncol(X)

  if(.Platform$OS.type == "windows" && parallel == TRUE){
    n.cores <- 1L
    parallel <- FALSE
    warning("Under Windows platform, parallel computing cannot be executed.")
  }

  #--------------------
  # Error Handling
  #--------------------
  if (!type %in% c("lasso", "grlasso", "slasso", "sgrlasso")) {
    stop("type has to be either lasso, grlasso, slasso or sgrlasso.")
  }

  if (!all(group==1:p) && (!type %in% c("grlasso", "sgrlasso"))) {
    stop("Choose type to be either grlasso or sgrlasso if group-structure exists.")
  }

  if (all(group==1:p) && (!type %in% c("lasso", "slasso"))) {
    stop("Choose type to be either lasso or slasso if group-structure does not exist.")
  }

  if (!PEtype %in% c("coeff", "mu")) {
    stop("Invalide PEtype input.")
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

  if (length(group) != p) {
    stop("group must have a same length with the number of X columns")
  }
  if (length(weights) != length(unique(group))) {
    stop("weights has to have a same length as the number of groups")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (any(!1:max(group) %in% group)) {
    stop("group index has to be a consecutive integer starting from 1.")
  }

  if (any(missing(pointEstimate_1), missing(sig2_1), missing(lbd_1))) {
    stop("provide all the parameters for the distribution")
  }
  # if (type == "lasso" && !sum(c(missing(sig2_2), missing(lbd_2),
  #                               missing(pointEstimate_2)))==3) {
  #   stop("Mixture distribution can be only used under group lasso.")
  # }

  if (!sum(c(missing(sig2_2), missing(lbd_2), missing(pointEstimate_2)))
      %in% c(0, 3)) {
    stop("provide all the parameters for the mixture distribution.")
  }
  #--------------------
  # Main Step
  #--------------------
  Mixture <- !missing(sig2_2)

  if (Mixture) {
    niter1 <- rbinom(n = 1, size = niter, prob = 1/2)
    niter2 <- niter-niter1

    PB1 <- PBsamplerMain(X = X, pointEstimate = pointEstimate_1, sig2 = sig2_1,
                         lbd = lbd_1, weights = weights, group = group, niter = niter1,
                         type = type, PEtype = PEtype, parallel = parallel, ncores = ncores, verbose = verbose)
    PB2 <- PBsamplerMain(X = X, pointEstimate = pointEstimate_2, sig2 = sig2_2,
                         lbd = lbd_2, weights = weights, group = group, niter = niter2,
                         type = type, PEtype = PEtype, parallel = parallel, ncores = ncores, verbose = verbose)
    if (type %in% c("lasso", "grlasso")) {
      RESULT <- list(beta = rbind(PB1$beta, PB2$beta),
                     subgrad = rbind(PB1$subgrad, PB2$subgrad), X = X,
                     pointEstimate = rbind(pointEstimate_1, pointEstimate_2),
                     sig2 = c(sig2_1, sig2_2), lbd = c(lbd_1, lbd_2), weights = weights, group = group,
                     type = type, PEtype = PEtype, mixture = Mixture)
    } else {
      RESULT <- list(beta = rbind(PB1$beta, PB2$beta),
                     subgrad = rbind(PB1$subgrad, PB2$subgrad), hsigma = c(PB1$hsigma, PB2$hsigma), X = X,
                     pointEstimate = rbind(pointEstimate_1, pointEstimate_2),
                     sig2 = c(sig2_1, sig2_2), lbd = c(lbd_1, lbd_2), weights = weights, group = group,
                     type = type, PEtype = PEtype, mixture = Mixture)
    }
  } else {
    PB <- PBsamplerMain(X = X, pointEstimate = pointEstimate_1,
      sig2 = sig2_1, lbd = lbd_1, weights = weights, group = group,
      niter = niter, type = type, PEtype = PEtype, parallel = parallel,
      ncores = ncores, verbose = verbose)
    if (type %in% c("lasso", "grlasso")) {
      RESULT <- list(beta = PB$beta, subgrad = PB$subgrad, X = X,
                     pointEstimate = pointEstimate_1, sig2 = sig2_1, lbd = lbd_1, weights = weights, group = group,
                     type = type, PEtype = PEtype, mixture = Mixture)
    } else {
      RESULT <- list(beta = PB$beta, subgrad = PB$subgrad, hsigma = PB$hsigma, X = X,
                     pointEstimate = pointEstimate_1, sig2 = sig2_1, lbd = lbd_1, weights = weights, group = group,
                     type = type, PEtype = PEtype, mixture = Mixture)
    }

  }
  class(RESULT) <- "PB"
  return(RESULT)
}

PBsamplerMain <- function(X, pointEstimate, sig2, lbd, weights = rep(1, max(group)),
 group, niter, type, PEtype, parallel,
 ncores, verbose)
{
  n <- nrow(X);
  p <- ncol(X);

  # Error Handling

  if (length(sig2) !=1 || length(lbd) != 1) {
    stop("sig2/lbd should be a scalar.")
  }

  if (sig2 <= 0) {
    stop("sig2 should be positive.")
  }

  if (lbd < 0) {
    stop("lbd should be non-negative.")
  }

  if (parallel && verbose) {
    warning("Note that verbose only works when parallel=FALSE")
  }

  if (PEtype == "coeff" && length(pointEstimate) != p) {
    stop("pointEstimate must have a same length with the col-number of X, if PEtype = \"coeff\"")
  }

  if (PEtype == "mu" && length(pointEstimate) != n) {
    stop("pointEstimate must have a same length with the row-number of X, if PEtype = \"mu\"")
  }

  if(verbose) {
    niter.seq <- round(seq(1, niter, length.out = 11))[-1]
  }

  W <- rep(weights, table(group))

  X.tilde   <- scale(X, FALSE, scale = W)
  t.X.tilde <- t(X.tilde) # same as  t(X) / W
  GramMat   <- t(X.tilde) %*% X.tilde
  Lassobeta <- matrix(0, niter, p)
  Subgrad   <- matrix(0, niter, p)
  Ysample   <- numeric(n)
  if (PEtype == "coeff") {
    Yexpect <- X %*% pointEstimate
  } else {
    Yexpect <- pointEstimate
  }
  if (type %in% c("lasso", "grlasso")) {
    FF <- function(x) {
      epsilon <- rnorm(n, mean = 0, sd = sig2^0.5)
      Ysample <- Yexpect + epsilon;
      #if(center){Ysample=Ysample-mean(Ysample);}

      LassoFit <- gglasso(X.tilde, Ysample, pf = rep(1, max(group)),
                          group = group, loss = "ls", intercept = FALSE, lambda = lbd)
      Lassobeta <- coef(LassoFit)[-1] / W
      return(c(Lassobeta, (t.X.tilde %*% Ysample -
                             GramMat %*% (Lassobeta * W)) / n / lbd))
    }
  } else {
    FF <- function(x) {
      epsilon <- rnorm(n, mean = 0, sd = sig2^0.5)
      Ysample <- Yexpect + epsilon;
      #if(center){Ysample=Ysample-mean(Ysample);}

      sig <- signew <- .5
      K <- 1 ; niter <- 0
      while(K == 1 & niter < 1000){
        sig <- signew;
        lam <- lbd * sig
        B0 <- coef(gglasso::gglasso(X.tilde,Ysample,loss="ls",group=group,pf=rep(1,max(group)),lambda=lam,intercept = FALSE))[-1]
        signew <- sqrt(crossprod(Ysample-X.tilde %*% B0) / n)

        niter <- niter + 1
        if (abs(signew - sig) < 1e-04) {K <- 0}
        if (verbose) {
          cat(niter, "\t", sprintf("%.3f", slassoLoss(X.tilde,Ysample,B0,signew,lbd)),"\t",
              sprintf("%.3f", signew), "\n")
        }
      }
      hsigma <- c(signew)
      S0 <- (t.X.tilde %*% (Ysample - X.tilde %*% B0)) / n / lbd / hsigma
      B0 <- B0 / rep(weights,table(group))
      return(c(B0, S0, hsigma))
    }
  }

  if (type %in% c("lasso", "grlasso")) {
    if (parallel == FALSE) {
      TEMP      <- lapply(1:niter, FF)
      TEMP      <- do.call(rbind, TEMP)
      Lassobeta <- TEMP[, 1:p]
      Subgrad   <- TEMP[, 1:p + p]
    } else {
      TEMP <- parallel::mclapply(1:niter, FF, mc.cores = ncores)
      TEMP <- do.call(rbind, TEMP)
      Lassobeta <- TEMP[, 1:p]
      Subgrad <- TEMP[, 1:p + p]
    }
    return(list(beta = TEMP[, 1:p], subgrad = TEMP[, 1:p + p]));
  } else {
    if (parallel == FALSE) {
      TEMP      <- lapply(1:niter, FF)
      TEMP      <- do.call(rbind, TEMP)
    } else {
      TEMP <- parallel::mclapply(1:niter, FF, mc.cores = ncores)
      TEMP <- do.call(rbind, TEMP)
    }
    return(list(beta = TEMP[, 1:p], subgrad = TEMP[, 1:p + p], hsigma = TEMP[, 2*p + 1]));
  }
}

#' @title Provide \code{(1-alpha)%} confidence interval of each coefficients
#'
#' @description Using samples drawn by \code{\link{PBsampler}}, compute
#' \code{(1-alpha)%} confidence interval of each coefficients.
#'
#' @param PBsample Bootstrap samples of class \code{PB} from \code{PBsampler}
#' @param alpha significance level.
#' @param method bias-correction method. Either to be "none" or "debias".
#' @param parallel Logical. If \code{TRUE}, use parallelization.
#' @param ncores Integer. The number of cores to use for the parallelization.
#'
#' @details If \code{method==none}, \code{\link{PB.CI}} simply compute
#' the quantile of the sampled coefficients. If \code{method==debias}, we use
#' debiased estimator to compute confidence interval.
#'
#' @return \code{(1-alpha)%} confidence interval of each coefficients
#'
#' @references
#' Zhang, C., Zhang, S. (2014) Confidence intervals for low dimensional
#' parameters in high dimensional linear models. Journal of the Royal
#' Statistical Society: Series B 76, 217–242.
#'
#' @examples
set.seed(1234)
n <- 40
p <- 50
Niter <-  10
Group <- rep(1:(p/10), each = 10)
Weights <- rep(1, p/10)
X <- matrix(rnorm(n*p), n)
object <- PBsampler(X,c(1,1,rep(0,p-2)),1,.5,niter=1000,type="lasso", parallel=TRUE, ncores = 8)
PB.CI(object = object, alpha = .05, method = "debias", parallel = TRUE, ncores = 8)

#' @export
PB.CI <- function(object, alpha = .05, method = "debias", parallel=TRUE, ncores=2L) {

  if (class(object)!="PB") {
    stop("object class has to be \"PB\".")
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

  X <- object$X
  weights <- object$weights
  lbd <- object$lbd
  group <- object$group

  if (!method %in% c("none", "debias")) {
    stop("method should be either \"default\" or \"debias\".")
  }
  if (method == "none") {
    Result <- apply(object$beta,2,quantile,prob=c(alpha/2,1-alpha/2))
  } else { # debiased estimator
    B <- object$beta # niter x p matrix
    S <- object$subgrad # niter x p matrix
    refitB <- matrix(0,nrow(B),ncol(B))
    W <- rep(weights, table(group))

    refitY <- solve(X%*%t(X))%*%X %*% (crossprod(X) %*% t(B) / n +
                                       lbd * W * t(S)) * n # n x niter matrix

    # Compute Z matrix, see Zhang, C. H. and Zhang S. S. (2014) JRSSB
    nodewiselasso.out <- score.nodewiselasso(x = X,
                                             parallel = parallel,
                                             ncores = ncores)
    Z <- nodewiselasso.out$out$Z
    scaleZ <- nodewiselasso.out$out$scaleZ
    Z <- scale(Z,center=FALSE,scale=1/scaleZ)
    #hdiFit <- hdi::lasso.proj(x = X, y = refitY[,1], standardize = FALSE, parallel = TRUE, return.Z = TRUE)
    #Z <- hdiFit$Z
    #hdiFit$bhat

    ZY <- crossprod(Z, refitY) # p x niter
    ZX <- colSums(Z * X) # length p
    ZXcomp <- matrix(0,p-1,p)
    for (i in 1:p) {
      ZXcomp[,i] <- crossprod(Z[,i], X[, -i])
    }

    refitB <- matrix(0,nrow(B),ncol(B))

    FF <- function(x) {
      TEMP <- c()
      for (j in 1:ncol(B)) {
        TEMP[j] <- (ZY[j,x] - ZXcomp[,j] %*% B[x,-j]) / ZX[j]
      }
      return(TEMP)
    }

    if (parallel) {
      refitB <- mcmapply(FF, 1:nrow(B), mc.cores = ncores)
    } else {
      refitB <- mapply(FF, 1:nrow(B))
    }
    # for (i in 1:nrow(B)) {
    #   for (j in 1:ncol(B)) {
    #     refitB[i,j] <- (ZY[j,i] - ZXcomp[,j] %*% B[i,-j]) / ZX[j]
    #   }
    # }
    Result <- apply(refitB,1,quantile,prob=c(alpha/2,1-alpha/2))
  }
  colnames(Result) <- paste("beta", 1:ncol(object$X), sep = "")
  #class(Result) <- "PB.CI"
  return(Result)
}
# HDI$bhat[j] == HDI$betahat[j] + t(Z[,j])%*%(Y-X%*%HDI$betahat) / crossprod(Z[,j],X[,j])
