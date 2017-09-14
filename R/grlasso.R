#' @title Compute group lasso estimator when lambda is given.
#'
#' @description Computes group lasso solution when lambda is given.
#'
#' @param X predictor matrix.
#' @param Y response vector.
#' @param lbd numeric. penalty term of group lasso.
#' @param weights weight vector with length equal to the number of groups. Default is
#' \code{weights = rep(1, max(group))}.
#' @param group \code{p} x \code{1} vector of consecutive integers describing the group structure.
#' The number of groups should be the same as max(group). Default is \code{group = 1:p}
#' , where \code{p} is number of covariates.
#' @param Gamma groupwise maximum eigenvalue multiplied by \code{(1 + 1e-05)}.
#' @param eps numeric. Tolerance level. Default is \code{eps = 1e-05}
#' @param returnGamma logical. If \code{returnGamma = TRUE}, return \code{Gamma}.
#' See \code{\link{cv.lasso}} for details.
#' @details
#' Computes group lasso solution.
#' When users use \code{\link{grlassoFit}} repetitively with same \code{X} and
#' \code{group}, they can accelerate the computational speed by providing \code{Gamma}.
#' \code{Gamma} can be easily obtained by setting \code{returmGamma = TRUE}.
#'
#' @references
#' Yang, Y. and Zou, H. (2015), “A Fast Unified Algorithm for Computing
#' Group-Lasso Penalized Learning Problems,” Statistics and Computing, 25(6), 1129-1141.
#'
#' @return \item{coef}{coefficient estimator.}
#' @return \item{Gamma}{see input arguments description.}
#' @examples
#' set.seed(1234)
#' n <- 10
#' p <- 20
#' group <- rep(1:4,each=5)
#'
#' X <- matrix(rnorm(n*p), n)
#' Y <- X %*% c(rep(1,5),rep(0, p-5)) + rnorm(n)
#' weights <- rep(1,max(group))
#' lbd <- .37
#'
#' grlassoFit(X = X, Y = Y, lbd = lbd, weights= weights, group = group,
#'  returnGamma = TRUE)
#'
#' Gamma <- grlassoFit(X = X, Y = Y, lbd = lbd, weights= weights, group = group,
#'  returnGamma = TRUE)$Gamma
#'
#' grlassoFit(X = X, Y = Y, lbd = lbd, weights= weights, group = group,
#'            Gamma = Gamma)
#' @export
grlassoFit <- function(X, Y, lbd, weights = rep(1, max(group)), group = 1:p,
                       Gamma, eps=1e-5, returnGamma = FALSE)
{
  n <- nrow(X)
  p <- ncol(X)
  #--------------------
  # Error Handling
  #--------------------
  if (lbd < 0) {
    stop("lbd should be non-negative.")
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
  if (!is.null(Y) & length(Y) != n) {
    stop("dimension of X and Y are not conformable.")
  }
  #--------------------
  XY <- X * c(Y)

  if (missing(Gamma)) {
    Gamma <- groupMaxEigen(X, group)
  }

  if (returnGamma) {
    return(list(coef = c(grlasso(X = X, Y = Y, XY = XY, weights = weights, group = group,
                                 lbd = lbd, Gamma = Gamma, initBeta = rep(1, ncol(X)),
                                 eps = eps))
    , Gamma = Gamma))
  } else {
    return(list(coef = c(grlasso(X = X, Y = Y, XY = XY, weights = weights, group = group,
                                 lbd = lbd, Gamma = Gamma, initBeta = rep(1, ncol(X)),
                                 eps = eps))))
  }
}
