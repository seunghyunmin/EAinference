#' @title MH sampler for lasso / group lasso estimator
#'
#' Metropolis-Hastings sampler for lasso / group lasso estimator
#' using estimator augmentation.
#'
#' @param X \code{n} x \code{p} matrix of predictors, where \code{n} is the
#' number of samples and \code{p} is the number of covariates.
#' @param pointEstimate numeric vector. Estimates of true coefficient for
#' \code{type="coeff"} or E(y) for \code{type="mu"}.
#' @param sig2 numeric. variance of error term.
#' @param lbd numeric. penalty term of lasso. See the loss function given at details.
#' @param group \code{p} x \code{1} vector. Consecutive integers should be used for indexing groups.
#' The number of groups should be same as \code{max(group)}. See examples for details.
#' @param weights \code{p} x \code{1} vector. weight term for each group.
#' @param B0 \code{p} x \code{1} vector. Initial value of lasso/group lasso estimator.
#' @param S0 \code{p} x \code{1} vector. Initial value of subgradients. If not given, it will be generated in defauly way.
#' @param A numeric vector. Active group index. \code{which(B0 != 0)} has to be a subset of \code{A}.
#' @param tau \code{|A|} x \code{1} numeric vector. Variance parameter for proposal
#' distribution of active coefficients.
#' @param niter numeric. The number of iterations.
#' @param burnin numeric. The length of burin-in periods
#' @param updateS.itv numeric. Update subgradients in every \code{updateS.itv} iterations. Set this value larger than \code{niter} if one wants to skip updating subgradients.
#' @param type either to be "coeff" or "mu". Decide what kind of \code{pointEstimate} to use.
#' @param verbose verbose
#' @param ... complementary arguments for MH-sampler for lasso.
#' \itemize{
#'  \item{\code{FlipSA}}{the parameter that is needed for the high-dimensional setting.
#' Has to be a subset of active set, A. If the index is not listed in FlipSA,
#' the sign of coefficients which corresponds to the index will be fixed.
#' The default is \code{FlipSA=A}}
#'  \item{\code{SFindex} }{subgradient index for the free coordinate.}
#'  \item{\code{randomSFindex} }{logical. If \code{true}, resample \code{SFindex} in every
#' \code{updateSF.itv} number.}
#'  \item{\code{updateSF.itv} }{Specifies how many iterations will be done without
#' updating the \code{SFindex}.}
#' }
#'
#' @details If futype="normal", it generate
#' @return \code{MHLS} returns an object of \code{\link{class}} \code{c("MHLS", "Lasso")} or \code{c("MHLS", "GroupLasso")}
#' The functions \code{summary} and \code{plot} are used for a breif summary and generating plots.
#' \item{beta}{lasso / group lasso samples}
#' \item{subgrad}{subgradient samples}
#' \item{acceptHistory}{number of acceptance and proposed}
#' \item{niteration}{number of iteration}
#' \item{burnin}{length of burn-in period}
#' \item{pointEstimate, group, type}{same as function arguments}
#' @examples
#' #-------------------------
#' # Low dim
#' #-------------------------
#' set.seed(123)
#' n <- 10
#' p <- 5
#' X <- matrix(rnorm(n * p), n)
#' Y <- X %*% rep(1, p) + rnorm(n)
#' sigma2 <- 1
#' lbd <- .37
#' weights <- rep(1, p)
#' LassoResult <- Lasso.MHLS(X = X, Y = Y, lbd = lbd, weights = weights)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' MHLS(X = X, pointEstimate = rep(0, p), sig2 = 1, lbd = 1, group = 1:p,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "coeff")
#' MHLS(X = X, pointEstimate = rep(0, n), sig2 = 1, lbd = 1, group = 1:p,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "mu")
#'
#' Group <- c(1, 1, 2, 2, 2)
#' weights <- rep(1,max(Group))
#' LassoResult <- Lasso.MHLS(X = X, Y = Y, lbd = lbd, weights = weights,
#'                           group = Group)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' MHLS(X = X, pointEstimate = rep(0, p), sig2 = 1, lbd = 1, group = Group,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "coeff")
#' MHLS(X = X, pointEstimate = rep(0, n), sig2 = 1, lbd = 1, group = Group,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "mu")
#'
#' #-------------------------
#' # High dim
#' #-------------------------
#' set.seed(123)
#' n <- 5
#' p <- 10
#' X <- matrix(rnorm(n*p),n)
#' Y <- X %*% rep(1,p) + rnorm(n)
#' sigma2 <- 1
#' lbd <- .37
#' weights <- rep(1,p)
#' LassoResult <- Lasso.MHLS(X = X,Y = Y,lbd = lbd,weights = weights)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' MHLS(X = X, pointEstimate = rep(0, p), sig2 = 1, lbd = 1, group = 1:p,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "coeff")
#' MHLS(X = X, pointEstimate = rep(0, n), sig2 = 1, lbd = 1, group = 1:p,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "mu")
#' @export
MHLS <-  function(X, pointEstimate, sig2, lbd, group = 1:ncol(X),
                   weights = rep(1, ncol(X)), B0, S0, A = unique(group[which(B0 != 0)]),
                   tau = rep(1, length(A)), niter = 2000, burnin = 0, type = "coeff", updateS.itv = 1, verbose = FALSE, ...)
{
  MHLSmain(X = X, pointEstimate = pointEstimate, sig2 = sig2, lbd = lbd,
    group = group, weights = weights, B0 = B0, S0 = S0, A = A,
           tau = tau, niter = niter, burnin = burnin, type = type, updateS.itv = updateS.itv, verbose = verbose, ...)
}

MHLSmain <- function (X, pointEstimate, sig2, lbd, group,
  weights, B0, S0, A, tau, niter, burnin, type, updateS.itv, verbose, ...)
{
  #------------------
  # Error handling
  #------------------
  if (type == "coeff" && length(pointEstimate) != p) {
    stop("pointEstimate must have a same length with the col-number of X, if type = \"coeff\"")
  }
  if (type == "mu" && length(pointEstimate) != n) {
    stop("pointEstimate must have a same length with the row-number of X, if type = \"mu\"")
  }
  if (!type %in% c("coeff", "mu")) {
    stop("Invalide type input.")
  }
  if (length(group) != p) {
    stop("group must have a same length with the number of X columns")
  }
  if (length(weights) != length(unique(group))) {
    stop("weights has to have a same length as the number of groups")
  }
  if (any(weights < 0)) {
    stop("weights should be non-negative.")
  }
  if (any(!1:max(group) %in% group)) {
    stop("group index has to be a consecutive integer starting from 1.")
  }
  if (any(missing(pointEstimate), missing(sig2), missing(lbd))) {
    stop("provide all the parameters for the distribution")
  }
  if (burnin >= niter) {
    stop("burnin has to be greater than niter")
  }
  if (niter <= 1) {
    stop("niter should be a integer greater than 1.")
  }
  if (all(group == 1:ncol(X))) {
    est <- MHLSswp(X = X, pointEstimate = pointEstimate, sig2 = sig2,
      lbd = lbd, weights = weights, B0 = B0, S0 = S0, A = A,
      tau = tau, niter = niter, burnin = burnin, type = type, updateS.itv, ...)
    class(est) <- "MHLS"
    class(est) <- append(class(est), "Lasso")
  } else {
    if (n < p) {
      stop("Low dimensional data can be only used for group lasso MH-sampler.")
    }
    est <- MHLSgroup(X = X, pointEstimate = pointEstimate, sig2 = sig2,
      lbd = lbd, weights = weights, group = group, B0 = B0, S0 = S0, A = A,
      tau = tau, niter = niter, burnin = burnin, type = type, updateS.itv = updateS.itv, ...)
    class(est) <- "MHLS"
    class(est) <- append(class(est), "GroupLasso")

    Group.matrix <- matrix(0, max(group), p)
    Group.matrix[cbind(group,1:p)] <- 1
    est$beta <- (est$group.l2.norm %*% Group.matrix) * est$subgrad
    est[1:4] <- est[c(1,4,2,3)]
    names(est)[1:4] <- names(est)[c(1,4,2,3)]
    #Need to be updated!!!!!!!!!!!!!
  }
  #est$acceptHistory <- rbind(est$acceptHistory,est$acceptHistory[1,]/est$acceptHistory[2,])
  est$niteration <- niter
  est$burnin <- burnin
  est$pointEstimate <- pointEstimate
  est$group <- group
  est$type <- type
  return(est)
}

MHLSswp <- function(X, pointEstimate, sig2, lbd, weights = rep(1, ncol(X)),
  B0, S0, A = which(B0 != 0), tau = rep(1, length(A)), niter = 2000,
  burnin = 0, type = "coeff", FlipSA = A, SFindex,
  randomSFindex = TRUE, updateSF.itv = round(niter/20), updateS.itv = 1,
  verbose = FALSE, ...)
{
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  NN <- floor((1:10) * niter / 10) # verbose index

  A <- unique(A) # active set
  Ac <- setdiff(1:p,A) # inactive set
  nA <- length(A) # size of the active set
  nI <- length(Ac) # size of the inactive set
  C <- crossprod(X) / n #Gram matrix

  if (!all(which(B0 != 0) %in% A)) {
    stop("The active set, A, has to include every index of nonzero B0.")
  }
  if (any(!A %in% 1:p)) {
    stop("The active set, A, has to be a subset of 1:ncol(X)")
  }
  if (!missing(S0) && !all(round(S0[which(B0 != 0)], 3) == sign(B0[B0 != 0]))) {
    stop("Invalid S0. Leave S0 blank, if S0 is unknown.")
  }
  if (length(tau) != length(A)) {
    stop("tau must have a same length with the active set, A.")
  }
  if (n >= p) {   # Low-dim MH
    #precalculation
    #for(j in 1:p){X[,j]=X[,j]-mean(X[,j])}
    #If X to be centered, we need to recompute B0 and S0 using centered X.
    Cinv <- solve(C) #Inverse Gram matrix
    logdiagC <- log(diag(C))
    Vinv <- n / (2 * sig2) * Cinv # non exponential part of pdf of U
    if (type == "coeff") {
      CB <- C %*% pointEstimate    # pointEstimate : True beta
    } else {
      CB <- crossprod(X, pointEstimate) / n    # pointEstimate : True beta
    }

    lbdwgt <- lbd * weights
    loglbd <- log(lbdwgt)

    #initialization
    B <- matrix(0, niter, p) #Store Coefficient
    S <- matrix(0, niter, p) #Store Sub-grad
    B[1,] <- B0 #B0 initial value of beta

    if (missing(S0)) {
      S[1, A] <- sign(B[1, A])
      S[1, Ac] <- runif(length(Ac), -1, 1)
    } else {
      S[1, ] <- S0
    }

    Ucur <- C %*% B[1, ] + lbdwgt * S[1, ] - CB #Current U
    ldetDRatio <- 0
    if (nA >= 1) {negCAAinv <- -solve(C[A, A])} else {negCAAinv <- NULL}

    nAccept <- numeric(2)
    nProp <- (niter - burnin) * c(nA, nI)

    if (length(setdiff(FlipSA, A))!=0)
      stop("FlipSA has to be a subset of active set, A.")
    A2 <- setdiff(A,FlipSA)
    if (any(B0[A2]==0))
      stop("To fix the sign of beta_j, use non-zero B0_j.")

    if (length(A2) != 0) {
      LUbounds <- matrix(0, p, 2);
      LUbounds[B0 < 0, 1] <- -Inf;
      LUbounds[B0 > 0, 2] <- Inf;
    }


    for(t in 2:niter)
    {
      if(nA >= 1){
        if (length(FlipSA)!=0) {
          for (j in FlipSA) {
            b_prop <- rnorm(1, mean = B[t - 1, j], sd = tau[which(A == j)])
            s_prop <- sign(b_prop)
            DiffU <- (b_prop - B[t - 1, j]) * C[, j]
            DiffU[j] <- DiffU[j] + lbdwgt[j] * (s_prop - S[t - 1, j])
            Uprop <- Ucur + DiffU
            logMH <- -t(Ucur + Uprop) %*% Vinv %*% DiffU
            u <- runif(1)
            if(log(u) < logMH)
            {
              B[t, j] <- b_prop
              S[t, j] <- s_prop
              Ucur <- Uprop
              if (t > burnin) {nAccept[1] <- nAccept[1] + 1}
              #nAccept[2]=nAccept[2]+1
            }else{
              B[t, j] <- B[t - 1, j]
              S[t, j] <- S[t - 1, j]
            }
          }
        }

        if (length(A2)!=0) {
          for (j in A2) {
            b_prop <- msm::rtnorm(1, mean = B[t - 1, j], sd = tau[which(A == j)],
                             lower = LUbounds[j, 1],
                             upper = LUbounds[j, 2])
            Ccur <- pnorm(0,mean=B[t-1, j],sd=tau[which(A == j)],lower.tail=(B[t-1,j]<0),log.p=FALSE);
            Cnew <- pnorm(0,mean=b_prop,sd=tau[which(A == j)],lower.tail=(b_prop<0),log.p=FALSE);
            lqratio=log(Ccur/Cnew);


            DiffU <- (b_prop - B[t - 1, j]) * C[, j]
            Uprop <- Ucur + DiffU
            logMH <- -t(Ucur + Uprop) %*% Vinv %*% DiffU + lqratio
            u <- runif(1)
            if(log(u) < logMH)
            {
              B[t, j] <- b_prop
              S[t, j] <- s_prop
              Ucur <- Uprop
              if (t > burnin) {nAccept[1] <- nAccept[1] + 1}
              #nAccept[2]=nAccept[2]+1
            }else{
              B[t, j] <- B[t - 1, j]
              S[t, j] <- S[t - 1, j]
            }
          }
        }

      }
      if(nI >= 1 && (t %% updateS.itv == 0)){
        for(j in Ac)
        {
          s_prop <- runif(1, -1, 1)
          diffu <- lbdwgt[j] * (s_prop - S[t - 1, j])
          Uprop <- Ucur
          Uprop[j] <- Ucur[j] + diffu
          logMH <- -t(Ucur + Uprop) %*% Vinv[, j] * diffu
          u <- runif(1)
          if(log(u) < logMH)
          {
            S[t,j] <- s_prop
            Ucur <- Uprop
            if (t > burnin) {nAccept[2] <- nAccept[2] + 1}
            #nAccept[3]=nAccept[3]+1
          } else {
            S[t, j] <- S[t - 1, j]
          }
        }
      } else {
        S[t, Ac] <- S[t - 1, Ac]
      }
      if (sum(t == NN)==1) {
        aa <- which(NN==t)
        cat(paste("Updating : ", aa * 10,"%" ,sep = ""), "\n")
      }
    }
    #nAccept=nAccept/c((niter-1)*selectsize,nProp)
    #nAccept=nAccept/nProp
  }
  if (n < p) {
    #precalculation---------------------------
    #for (j in 1:p) {X[,j]=X[,j]-mean(X[,j])}
    #If X to be centered, we need to recompute B0 and S0 using centered X.
    C <- t(X) %*% X / n
    egC <- eigen(C)
    V <- egC$vectors
    R <- 1:n
    N <- (n + 1):p
    InvVarR <- 1 / (egC$values[R] * sig2 / n) #inverse of (sig2*Lambda_i/n)
    VR <- matrix(V[,R], p, n)
    VRC <- t(VR) %*% C
    W <- diag(weights)
    LBD <- diag(egC$values[R])
    lbdVRW <- lbd * t(VR) %*% W
    if (type == "coeff") {
      VRCB <- t(VR) %*% C %*% pointEstimate
    } else {
      VRCB <- t(VR) %*% crossprod(X, pointEstimate) / n
    }


    # if (is.missing(B0)) {
    #   if (signBA == NULL)
    #     stop("If initial value of 'B0' is not given, 'signBA' has to be given.")
    #   B0 = rep(0,p)
    #   B0[A] = abs(rnorm(nA)) * signBA
    # }

    if (length(setdiff(FlipSA, A))!=0)
      stop("FlipSA has to be a subset of active set, A.")
    A2 <- setdiff(A,FlipSA)

    V_AN <- V[A,N,drop=FALSE]
    V_IN <- V[Ac,N,drop=FALSE]
    V_AR <- V[A,R,drop=FALSE]
    V_IR <- V[Ac,R,drop=FALSE]

    BB <- t(V_IN) %*% W[Ac,Ac]
    tVAN.WA <- t(V_AN)%*%W[A,A]

    if (!missing(S0) && !all.equal(t(V[,N])%*%S0,matrix(0,length(N),1))) {
      warning("Invalid S0. Regenerate S0 with a default way")
      S0 <- NULL
    }

    if (missing(S0) || is.null(S0)) {
      E <- BB
      H <- rep(-1,2*nI)
      G <- rbind(diag(rep(1,nI)),diag(rep(-1,nI)))
      F1 <- -tVAN.WA%*%sign(B0[A])
      S0.prop <- limSolve::lsei(G=G,H=H,E=E,F=F1)
      if (S0.prop$IsError) {
        stop("There exist no solution for the given 'B0'.")
      } else {
        S0 <- rep(0,p)
        S0[A] <- sign(B0[A])
        S0[-A] <- S0.prop$X
      }
    }

    if (n-nA != 0) {
      # SFindex : S_F, free coordinate, index
      if (missing(SFindex)) {SFindex <- 1:(n-nA)}

      if (length(SFindex) != n-nA) {
        warning("Length of SFindex has to be same as 'n-length(A)'. Automatically set 'SFindex <- 1:(n-nA)'.")
        SFindex <- 1:(n-nA)
      }

      B_F <- BB[, SFindex, drop=FALSE]
      B_D <- BB[,-SFindex, drop=FALSE]
      invB_D <- solve(B_D)
      BDF <- invB_D%*%B_F
    } else if (n-nA == 0) {
      if (missing(SFindex)) {SFindex <- NULL}
      if (!is.null(SFindex)) {
        warning("If 'n-nA == 0', SFindex has to be set to NULL. Automatically set 'SFindex <- NULL'.")
        SFindex <- NULL
      }
      BDF <- invB_D <- solve(BB)  # NEED TO CHECK
      B_F <- NULL
    }

    #initialization--------------------------
    B <- matrix(0,niter,p)
    S <- matrix(0,niter,p)
    B[1,] <- B0
    S[1,] <- S0

    Rcur <- VRC[ , A, drop=FALSE] %*% t(B[1, A, drop=FALSE]) + lbdVRW %*% S[1, ] - VRCB
    logfRcur <- -0.5 * sum(Rcur^2 * InvVarR)
    Tcur <- CalTmat(p,n,V,LBD,W,lbd,R,N,A,Ac)
    logdensitycur <- logfRcur+Tcur$logdetT

    #record acceptance rates
    nAccept <- numeric(2)
    #nProp=numeric(2)
    nProp <- c(nA*(niter-max(1,burnin)),(n-nA)*(niter-max(1,burnin)))
    #Change sign count
    nSignChange <- numeric(3)
    for(t in 2:niter )
    {
      #P1: update b_A
      if(nA>=1)
      {
        S[t,] <- S[t-1,]
        if (length(FlipSA)!=0) {
          MoveBA <- UpdateBA(B[t-1,],S[t-1,],tau,A,Ac,Rcur,logfRcur,VRC,lbdVRW,InvVarR,
                          tVAN.WA,invB_D,B_F,FlipSA,SFindex)
          B[t,] <- MoveBA$B
          S[t,] <- MoveBA$S
          nSignChange <- nSignChange+MoveBA$nChangeSign
          Rcur <- MoveBA$Rvec
          logfRcur <- MoveBA$logf
          nAccept[1] <- nAccept[1]+MoveBA$nAccept
          #nProp[1]=nProp[1]+nA
        }

        if (length(A2)!=0) {
          MoveBA <- UpdateBA.fixedSA(B[t-1,],tau,A2,Rcur,logfRcur,VRC,InvVarR)
          B[t,A2] <- MoveBA$B[A2]
          #S[t,]=S[t-1,]
          Rcur <- MoveBA$Rvec
          logfRcur <- MoveBA$logf
          nAccept[1] <- nAccept[1] + MoveBA$nAccept
        }
      }else{ B[t,] <- B[t-1,]; S[t,] <- S[t-1,]}

      # P2: update S_I
      if(nA<n && (t %% updateS.itv == 0))
      {
        MoveSI <- UpdateSI(S[t,],A,Ac,Rcur,n,p,logfRcur,lbdVRW,InvVarR,
                        tVAN.WA,invB_D,BDF,B_F,SFindex,...)
        S[t,] <- MoveSI$S
        Rcur <- MoveSI$Rvec
        logfRcur <- MoveSI$logf
        nAccept[2] <- nAccept[2]+MoveSI$nAccept
        #nProp[2]=nProp[2]+(n-nA)
      }#else{S[t,]=S[t-1,];S[t,A]=sign(B[t,A])}

      if (!is.null(SFindex) && randomSFindex && (t%%updateSF.itv==0) ) {
        SFindex <- sort(sample(1:(p-nA),n-nA))
        if (verbose) {cat("New SFindex : [",paste(SFindex,collapse=", "),"]\n")}
        B_F <- BB[,SFindex,drop=FALSE]
        B_D <- BB[,-SFindex,drop=FALSE]
        invB_D <- solve(B_D)
        BDF <- invB_D%*%B_F
      }

      if (sum(t == NN)==1) {
        aa <- which(NN==t)
        cat(paste("Updating : ", aa*10  ,"%",sep=""),"\n")
      }
    }
  }
  return(list(beta = B[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
              subgrad = S[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
              acceptHistory = rbind(nAccept, nProp)))
  #
  # return(list(beta = B[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
  #             subgrad = S[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
  #             pointEstimate = pointEstimate,
  #             signchange=nSignChange,
  #             acceptHistory = rbind(nAccept, nProp)))
}

MHLSgroup <- function(X, pointEstimate, sig2, lbd,
 weights, group, B0, S0, A, tau, niter, burnin, type = "coeff", updateS.itv, verbose)
{
  K <- 10
  W <- rep(weights,table(group))
  Psi <- 1/n * crossprod(X)
  if (n > p) {
    inv.Var <- n/sig2 * solve(Psi)
  }
  nA <- length(A)

  r.seq <- matrix(, niter, max(group))
  S.seq <- matrix(, niter ,p)
  nAccept <- numeric(2)
  nProp <- c(nA*(niter-max(1,burnin)), max(group)*(niter-max(1,burnin)))

  rcur <- group.norm2(B0,group)
  r.seq[1, ] <- rcur
  S.seq[1, ] <- Scur <- S0

  if (type == "coeff") {
    Hcur <- drop(Psi %*% drop(B0 - pointEstimate) + lbd * W * drop(S0))
  } else {
    Hcur <- drop(Psi %*% drop(B0) - t(X) %*% pointEstimate / n + lbd * W * drop(S0))
  }

  if (n >= p) {
    for (i in 2:niter) {
      r.new <- ld.Update.r(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,tau,type)
      r.seq[i,] <- rcur <- r.new$r
      Hcur <- r.new$Hcur
      if (i > burnin) {nAccept[1] <- nAccept[1] + r.new$nrUpdate}

      if (i %% updateS.itv == 0) {
        S.new <- ld.Update.S (rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,p,type)
        S.seq[i,] <- Scur <- S.new$S
        Hcur <- S.new$Hcur
      } else {
        S.seq[i,] <- Scur
      }
      if (i > burnin) {nAccept[2] <- nAccept[2] + S.new$nSUpdate}

      if (verbose && (i %% round(niter/10) == 0)) {
        cat("MCMC step,", K, "% Finished\n")
        K <- K+10
      }
    }
  # } else {
  #   for (i in 2:niter) {
  #     r.new <- hd.Update.r(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,1)
  #     r.seq[i,] <- rcur <- r.new$r
  #     Hcur <- r.new$Hcur
  #     nAccept[1] <- nAccept[1] + r.new$nrUpdate
  #
  #     S.new <- hd.Update.S (rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,p)
  #     S.seq[i,] <- Scur <- S.new$S
  #     Hcur <- S.new$Hcur
  #     nAccept[2] <- nAccept[2] + S.new$nSUpdate
  #
  #     if (verbose && (i %% round(niter/10) == 0)) {
  #       cat("MCMC step,", K, "% Finished\n")
  #       K <- K+10
  #     }
  #   }
  }
  return(list(
    group.l2.norm = r.seq[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
    subgrad = S.seq[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
    acceptHistory = rbind(nAccept, nProp)))
}

#' @export
print.MHLS <- function (x) {
  cat ("===========================\n")
  cat ("Number of iteration: ", x$niteration,"\n\n")
  cat ("Burn-in period: ", x$burnin,"\n\n")
  cat ("Plug-in pointEstimate: \n")
  print(x$pointEstimate)
  cat ("type: \n")
  print(x$type)

  # if (inherits(x,"Group")) {
  #   Group.matrix <- matrix(0, length(unique(x$group)), p)
  #   Group.matrix[cbind(x$group,1:p)] <- 1
  #   beta <- (x$group.l2.norm %*% Group.matrix) * x$subgrad
  # }
  cat ("\nLast 10 steps of beta's:\n")

  if (inherits(x,"GroupLasso")) {
    if (x$niteration-x$burnin <= 9) {
      print(x$group.l2.norm)
    } else {
      print(x$group.l2.norm[(x$niteration-x$burnin-9):(x$niteration-x$burnin),])
    }
  }

  if (x$niteration-x$burnin <= 9) {
    # if (inherits(x,"Group")) {
    #   #print(x$group.l2.norm)
    #   print(beta)
    # } else {
    print(x$beta)
    # }
  } else {
    # if (inherits(x,"Group")) {
    #   print(beta[(x$niteration-x$burnin-9):(x$niteration-x$burnin),])
    #   #print(x$group.l2.norm[(x$niteration-x$burnin-9):(x$niteration-x$burnin),])
    # } else {
    print(x$beta[(x$niteration-x$burnin-9):(x$niteration-x$burnin),])
    # }
  }

  cat ("\nlast 10 steps of subgradients:\n")
  if (x$niteration-x$burnin <= 9) {
    print(x$subgrad)
  } else {
    print(x$subgrad[(x$niteration-x$burnin-9):(x$niteration-x$burnin),])
  }

  cat ("\nAcceptance rate:\n")
  cat("-----------------------------\n")
  if (inherits(x,"GroupLasso")) {
    cat("\t \t l_2 group norm\t subgrad\n")
  } else {
    cat("\t \t \t beta \t subgrad\n")
  }
  cat("# Accepted\t : \t", paste(x$acceptHistory[1,],"\t"),"\n")
  cat("# Moved\t\t : \t", paste(x$acceptHistory[2,],"\t"),"\n")
  cat("Acceptance rate\t : \t", paste(round(x$acceptHistory[1,]/x$acceptHistory[2,],3),"\t"),"\n")
  # cat ("\nSignChange rate:\n")
  # cat("-----------------------------\n")
  # cat("# Accepted\t : \t", paste(x$signchange[1],"\t"),"\n")
  # cat("# Moved\t\t : \t", paste(x$signchange[2],"\t"),"\n")
  # cat("# Cdt Accept \t : \t", paste(x$signchange[3],"\t"),"\n")
  # cat("Acceptance rate\t : \t", paste(round(x$signchange[1]/x$signchange[2],3),"\t"),"\n")
}


#' @title Summarizing Metropolis-Hastings sampler outputs
#'
#' @description summary method for class "MHLS"
#'
#' @param object an object of class "MHLS", which is a result of \code{\link{lm}}.
#' @details
#' This function provides a summary of each sampled beta and subgradient.
#' @return mean, median, s.d., 2.5% quantile and 97.5% quantile for each sampled beta and subgradient.
#' @examples
#' #' set.seed(123)
#' n <- 10
#' p <- 5
#' X <- matrix(rnorm(n * p), n)
#' Y <- X %*% rep(1, p) + rnorm(n)
#' sigma2 <- 1
#' lbd <- .37
#' weights <- rep(1, p)
#' LassoResult <- Lasso.MHLS(X = X, Y = Y, lbd = lbd, weights = weights)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' summary(MHLS(X = X, pointEstimate = rep(0, p), sig2 = 1, lbd = 1, group = 1:p,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "coeff"))
#' @export
summary.MHLS <- function ( object ) {
  betaSummary <- t(apply(object$beta,2,SummBeta))
  #signsummary <- t(apply(object$beta,2,SummSign))
  subgradSummary <- t(apply(object$subgrad,2,SummBeta))
  result <- list(beta=betaSummary,subgradient=subgradSummary)
  class(result) <- "summary.MHLS"
  return(result)
}

#' @title Plotting Metropolis-Hastings sampler outputs
#'
#' @description For each index, this provides six plots;
#'  histogram, path plot and acf plot for beta and subgradient.
#'
#' @param object an object of class "MHLS", which is a result of \code{\link{lm}}.
#' @param A an index of covariates that one can plot with.
#' @param skipS logical. If \code{TRUE}, plot betas only.
#' @details
#' \code{plot.MHLS} provides summary plots of each sampled beta and subgradient.
#'  The first column provides histogram of beta and subgradient, while the second
#'  and the third column provide path plot and acf plot, respectively.
#' @return this function gives mean, median, s.d., 2.5% quantile and 97.5%
#'  quantile for each sampled beta and subgradient.
#' @examples
#' #' set.seed(123)
#' n <- 10
#' p <- 5
#' X <- matrix(rnorm(n * p), n)
#' Y <- X %*% rep(1, p) + rnorm(n)
#' sigma2 <- 1
#' lbd <- .37
#' weights <- rep(1, p)
#' LassoResult <- Lasso.MHLS(X = X, Y = Y, lbd = lbd, weights = weights)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' plot(MHLS(X = X, pointEstimate = rep(0, p), sig2 = 1, lbd = 1, group = 1:p,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "coeff"))
#' @export
plot.MHLS <- function(object, A=1:ncol(object$beta), skipS=FALSE, ... ) {
  #	n=nrow(object$beta)
  niter <- object$niteration
  burnin <- object$burnin

  if (!skipS) {par(mfrow = c(2,3))} else {par(mfrow = c(1,3))}

  if (!skipS)	{
    for (i in A) {
      hist(object$beta[,i],breaks=20,prob=T,xlab=paste("Beta_",i,sep=""),ylab="Density",main="")
      #ts.plot(object$beta[,i],xlab="Iterations",ylab="Samples")
      plot((burnin+1):niter,object$beta[,i],xlab="Iterations",ylab="Path",type="l")
      if ( sum(abs(diff(object$beta[,i]))) == 0 ) { plot( 0,type="n",axes=F,xlab="",ylab="")
        text(1,0,"Auto correlation plot \n not available",cex=1)} else {
          acf(object$beta[,i],xlab="Lag",main="")
        }
      hist(object$subgrad[,i],breaks=seq(-1-1/10,1,by=1/10)+1/20,prob=T,xlim=c(-1-1/20,1+1/20),xlab=paste("Subgradient_",i,sep=""),ylab="Density",main="")
      #ts.plot(object$subgrad[,i],xlab=Iterations,ylab="Samples")
      plot((burnin+1):niter,object$subgrad[,i],xlab="Iterations",ylab="Samples",type="l")
      if ( sum(abs(diff(object$subgrad[,i]))) == 0 ) { plot( 0,type="n",axes=F,xlab="",ylab="")
        text(1,0,"Auto correlation plot \n not available",cex=1)} else {
          acf(object$subgrad[,i],xlab="Lag",main="")
        }
      readline("Hit <Return> to see the next plot: ")
    }
  } else {
    for (i in A) {
      hist(object$beta[,i],breaks=20,prob=T,xlab=paste("Beta_",i,sep=""),ylab="Density",main="")
      #ts.plot(object$beta[,i],xlab="Iterations",ylab="Samples")
      plot((burnin+1):niter,object$beta[,i],xlab="Iterations",ylab="Path",type="l")
      if ( sum(abs(diff(object$beta[,i]))) == 0 ) { plot( 0,type="n",axes=F,xlab="",ylab="")
        text(1,0,"Auto correlation plot \n not available",cex=1)} else {
          acf(object$beta[,i],xlab="Lag",main="")
        }
      readline("Hit <Return> to see the next plot: ")
    }
  }
}



