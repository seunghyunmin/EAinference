#' This function directly draw samples with bootstrap using gglasso package.
#'
#' @param X n x p matrix of predictors.
#' @param coeff n x 1 vector of estimates of true coefficient.
#' @param sigma2 variance of error term. Not needed for futype="nonparametric".
#' @param lbd penalty term of lasso. See the loss function given at details.
#' @param weights Weight term for each group. Default is rep(1, max(group)).
#' @param group p x 1 vector of consecutive integers. The number of groups should be same as max(group).
#' @param niter The number of iterations.
#' @param futype Method for resampling errors. Either to be "normal" or "nonparametric". Default is "normal".
#' @param Y n x 1 vector of response. Not needed for futype="nonparametric".
#' @param center Whether to center Y values. Default is FALSE.
#' @param parallel Whether to parallelize the code. Default is FALSE.
#' @param ncores The number of cores to use for the parallelization. If missing, it uses maximum number of cores.
#' @param verbose Whether to show the process. Default is FALSE. Only works when parallel=FALSE.
#'
#' @details If futype="normal", it generate
#' @return \describe{
#'   \item{beta}{coefficient matrix of size N x p.}
#'   \item{subgrad}{subgradient matrix of size N x p.}
#'  }
#' @examples
#' n <- 10
#' p <- 30
#' niter <-  10
#' group <- rep(1:(p/10), each=10)
#' Weights <- rep(1,p/10)
#' x <- matrix(rnorm(n*p), n)
#' DirectSampler(X = x, coeff = rep(0,p), sigma2=1, lbd=.5, weights=Weights, group=group,N=niter, parallel=FALSE)
#' DirectSampler(X = x, coeff = rep(0,p), sigma2=1, lbd=.5, weights=Weights, group=group,N=niter, parallel=TRUE)
#' @export
DirectSampler=function(X,coeff,sigma2,lbd,weight=rep(1,max(group)),group=1:ncol(X),niter=2000,futype="normal",Y,center=FALSE, parallel=FALSE, ncores, verbose=FALSE)
{
  n <- nrow(X);
  p <- ncol(X);

  if (missing(Y) && futype=="nonparametric") {
    stop("Y is needed for the futype=\"nonparametric\"")
  }
  if (parallel && verbose) {
    warning("Note that verbose only works when parallel=FALSE")
  }
  if (length(coeff) != p) {
    stop("coeff must have a same length with the number of X columns")
  }
  if (length(group) != p) {
    stop("group must have a same length with the number of X columns")
  }
  if (length(weights) != length(unique(group))) {
    stop("weights has to have a same length as the number of groups")
  }
  if (futype %in% c("normal","nonparametric") == FALSE){
    stop("futype needs to be either normal or nonparametric")
  }
  if (any(!1:max(group) %in% group)) {
    stop("Group index has to be a consecutive integer.")
  }


  if(verbose) {
    ptm <- proc.time();
    niter.seq <- round(seq(1,niter,length.out = 11))[-1]
  }

  W <- rep(weights, table(group))

  X.tilde   <- scale(X,FALSE,scale=W)
  t.X.tilde <- t(X.tilde) # same as  t(X) / W
  GramMat   <- t(X.tilde) %*% X.tilde;
  Lassobeta <- matrix(0, niter, p);
  Subgrad   <- matrix(0, niter, p);
  Ysample   <- numeric(n);
  Yexpect   <- X %*% coeff;

  if(futype=="nonparametric"){obsRed=Y-Yexpect;}

  FF = function(x) {
    if(futype=="normal")
    {
      epsilon=rnorm(n,mean=0,sd=sigma2^0.5);
    }
    if(futype=="nonparametric")
    {
      epsilon=sample(obsRed,size=n,replace=TRUE);
    }
    Ysample=Yexpect+epsilon;
    if(center){Ysample=Ysample-mean(Ysample);}

    LassoFit <- gglasso(X.tilde, Ysample, pf = rep(1,max(group)), group = group, loss="ls", intercept=F, lambda=lbd)
    Lassobeta <- coef(LassoFit)[-1] / W
    return(c(Lassobeta,(t.X.tilde%*%Ysample-GramMat%*%(Lassobeta * W)) / n / lbd))
  }

  if (parallel == FALSE) {
    TEMP <- lapply(1:niter,FF)
    TEMP <- do.call(rbind,TEMP)

    Lassobeta <- TEMP[,1:p]
    Subgrad <- TEMP[,-c(1:p)]
  } else {
    if (missing(ncores)) {ncores <- parallel::detectCores()}
    options(mc.cores=ncores)

    TEMP <- parallel::mclapply(1:niter,FF)
    TEMP <- do.call(rbind,TEMP)

    Lassobeta <- TEMP[,1:p]
    Subgrad <- TEMP[,-c(1:p)]
  }
  if(verbose) {print(proc.time() - ptm)}
  return(list(beta=Lassobeta,subgrad=Subgrad));
}
