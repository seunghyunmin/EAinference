#' @title cross validation for scaled lasso & scaled group lasso
#'
#' @description Computes K-fold cross-validation to propose a lambda value
#'
#' @param X predictor matrix of size \code{p} x \code{n}.
#' @param Y response vector of length \code{n}.
#' @param group p x 1 vector of consecutive integers. The number of groups
#' should be same as max(group).
#' @param weights weight vector of length \code{max(group)}.
#' Weight term for each group. Default is \code{rep(1, max(group))}.
#' @param type type of penalty, either to be "lasso", "grlasso", "slasso" or "sgrlasso".
#' @param K integer. number of folds
#' @param minlbd numeric. Minumum value of the lambda sequence.
#' @param maxlbd numeric. Maximum value of the lambda sequence.
#' @param num.lbdseq integer. Number of the lambda sequence.
#' @param parallel logical. If \code{TRUE}, use parallelization.
#' @param ncores integer. The number of cores to use for the parallelization.
#' @param plot.it logical. If ture, plot the squared error curve
#' @param verbose verbose
#'
#' @details Use K-fold cross-validaiton to propose a optimal lambda.
#'
#' @return \item{lbd.min}{a value of lambda such gives a minimum squared error.}
#' @return \item{lbd.1se}{a largest lambda within 1 std. from \code{lbd.min}.}
#' @return \item{lbd.seq}{lambda sequence.}
#' @return \item{cv}{mean squared error at each lambda value.}
#' @return \item{cvsd}{the standard deviaion for cv.}
#'
#' @examples
#' set.seed(123)
#' n <- 30
#' p <- 50
#' group <- rep(1:(p/10),each=10)
#' weights <- rep(1, max(group))
#' X <- matrix(rnorm(n*p),n)
#' truebeta <- c(rep(1,5),rep(0,p-5))
#' Y <- X%*%truebeta + rnorm(n)
#'
#' cv.lasso(X,Y,group,weights,K=3,type="grslasso",num.lbdseq=10,plot=TRUE)
#' cv.lasso(X,Y,group,weights,K=10,type="grlasso",num.lbdseq=100,plot=TRUE)
#' @export
cv.lasso <- function(
  X,
  Y,
  group = 1:ncol(X),
  weights = rep(1,max(group)),
  type,
  K = 10L,
  minlbd = 0,
  maxlbd = if(type %in% c("lasso","grlasso")){2}else{max(t(X)%*%Y)},
  num.lbdseq = 300L,
  parallel = FALSE,
  ncores = 2L,
  plot.it = FALSE,
  verbose = FALSE)
{
  n <- nrow(X)
  p <- ncol(X)

  Y <- as.vector(Y)
  K <- as.integer(K)
  num.lbdseq <- as.integer(num.integer)

  #--------------------
  # Error Handling
  #--------------------
  if (n < 10) {stop("Sample size is too small to run cross-validation.")}
  if (n < 30) {warning("Insufficient sample size. The result can be unstable.")}

  if (length(Y) != n) {
    stop("dimension of X and Y are not conformable.")
  }
  if (!type %in% c("lasso", "grlasso", "slasso", "sgrlasso")) {
    stop("type has to be either lasso, grlasso, slasso or sgrlasso.")
  }
  if (length(group) != p) {
    stop("group must have a same length with the number of X columns")
  }
  if (!all(group==1:p) && (!type %in% c("grlasso", "sgrlasso"))) {
    stop("Choose type to be either grlasso or sgrlasso if group-structure exists.")
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
  if (length(weights) != length(unique(group))) {
    stop("weights has to have a same length as the number of groups")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (any(!1:max(group) %in% group)) {
    stop("group index has to be a consecutive integer starting from 1.")
  }
  if (any(c(minlbd,maxlbd) < 0)) {stop("minlbd/maxlbd should be non-negative")}
  if (num.lbdseq <= 0) {
    stop("num.lbdseq should be non-negative")
  }

  all.folds <- split(sample(1:n), rep(1:K, length = n))

  #index=seq(0,max(t(X)%*%Y),length=100)
  index <- seq(minlbd,maxlbd,length=num.lbdseq)[-1]
  residmat <- matrix(0, length(index), K)

  FF <- function(x,omit) {
    fit <- Lasso.MHLS(X[-omit,,drop=FALSE],Y[-omit],type=type,
                      lbd=index[x],group=group,weights=weights)$B0
    fit <- X[omit,,drop=FALSE]%*%fit
    return(mean((Y[omit]-fit)^2))
  }

  for (i in seq(K)) {
    omit <- all.folds[[i]]
    residmat[,i] <- do.call(rbind,parallel::mclapply(1:length(index),
      FF,omit = omit, mc.cores = ncores))
    if(verbose) {cat("\n CV Fold", i, "\n\n")}
  }
  #apply(residmat, 2, which.min)
  cv <- apply(residmat, 1, mean)
  cvsd <- apply(residmat, 1, sd)
  index.min.cv <- which.min(cv)

  err.1se <- cvsd[index.min.cv] + cv[index.min.cv]

  lbd.1se <- index[which.min(abs(cv - err.1se))]
  lbd.min <- index[index.min.cv]

  if (plot.it) {
    matplot(x=index, y=cbind(cv,cv-cvsd,cv+cvsd), type="l",
            lty=c(1,2,2), col=c(1,2,2),
            xlab="lambda",ylab="squared error",
            main="cross validation")
    abline(v=lbd.min,lty=2)
    abline(v=lbd.1se,lty=3)
  }
  return(list(lbd.seq=index, cv=cv, cvsd=cvsd,
              lbd.min=lbd.min, lbd.1se=lbd.1se,
              type=type))
}
