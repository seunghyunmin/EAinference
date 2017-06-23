#' Metropolis-Hastings sampler using estimator augmentation in lasso estimator.
#'
#' @param X n x p matrix of predictors.
#' @param pointEstimate n x 1 vector of estimates of true coefficient.
#' @param sig2 variance of error term.
#' @param weights Weight term for each coefficient. Default is rep(1, p).
#' @param lbd penalty term of lasso. See the loss function given at details.
#' @param group p x 1 vector of consecutive integers.
#' The number of groups should be same as max(group).
#' @param niter The number of iterations. Default value is 2000.
#' @param burnin The length of burin-in periods. Default value is 0.
#' @param B0 Initial value of coefficients of length p.
#' Make sure to set values to zero for inactive set and non-zero for active set.
#' @param S0 Initial value of subgradients. Has to satisfy certain conditions
#' for subgradients. If not given, it will be generated in defauly way
#' @param tau p x 1 vector. Variance parameter for proposal
#' distribution of coefficients.
#' @param indexweights p x 1 vector. Weights for updating active/inactive Index.
#' @param FlipSA The parameter that is needed for the high-dimensional setting.
#' Has to be a subset of active set, A. If the index is not listed in FlipSA,
#' the sign of coefficients which corresponds to the index will be fixed.
#' The default is \code{FlipSA=A}.
#' @param SFindex Subgradient index for the free coordinate.
#' @param randomSFindex logical. If \code{true}, resample \code{SFindex} in every
#' \code{updateSFindex} number.
#' @param updateSFindex Specifies how many iterations will be done without
#' updating the SFindex.
#' @param verbose verbose
#' @param updateS.itv Update S in every updateS.itv iterations. Set this value larger than niter if one wants to skip updating subgradients.
#'
#' @details If futype="normal", it generate
#' @return \describe{
#'   \item{beta}{coefficient matrix of size N x p.}
#'   \item{subgrad}{subgradient matrix of size N x p.}
#'  }
#' @examples
#' n <- 10
#' p <- 30
#' lbd <- .5
#' niter <-  10
#' group <- rep(1:(p/10), each=10)
#' Weights <- rep(1,p/10)
#' x <- matrix(rnorm(n*p), n)
#' DirectSampler(X = x, pointEstimate = rep(0,p), sig2=1, lbd=.5, weights=Weights,
#' group=group,N=niter, parallel=FALSE)
#' DirectSampler(X = x, pointEstimate = rep(0,p), sig2=1, lbd=.5, weights=Weights,
#' group=group,N=niter, parallel=TRUE)
#' @export
MHLS <- function (X, pointEstimate, sig2, lbd, group = 1:ncol(X),
                  weights = rep(1, ncol(X)), B0, S0, A = which(B0!=0), tau = rep(1, length(A)),
                  niter=2000, burnin=0, type = "coeff", ...)
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
    est <- MHLSswp(X,pointEstimate,sig2,lbd,weights,B0,S0,A,tau,niter,burnin,type,
                    ...)
  } else {
    est <- fixedA.MCMC(X,pointEstimate,sig2,weights,lbd,group,niter=2000,A,B0,S0,tau,verbose=FALSE)
    #Need to be updated!!!!!!!!!!!!!
  }
  est$niteration <- niter
  est$burnin <- burnin
  est$pluginbeta <- pointEstimate
  est$type <- type
  class(est) <- "MHLS"
  return(est)
}


#MHLS <-  function (X , ...) UseMethod("MHLS")

MHLSswp <- function(X, pointEstimate, sig2, lbd, weights = rep(1, ncol(X)),
  B0, S0, A = which(B0 != 0), tau = rep(1, length(A)), niter = 2000,
  burnin = 0, type = "coeff", FlipSA = A, SFindex,
  randomSFindex = TRUE, updateSFindex = round(niter/20), updateS.itv = 1,
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

  if (type == "coeff") {
    if (any(c(length(pointEstimate),length(weights)) != p)) {
      stop("pointEstimate/weights must have a same length with the number of X columns.")
    }
  } else {
    if (any(c(length(pointEstimate),length(weights)) != p)) {
      stop("pointEstimate/weights must have a same length with the number of X columns.")
    }
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
    CB <- C %*% pointEstimate    # pointEstimate : True beta
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
            b_prop <- rtnorm(1, mean = B[t - 1, j], sd = tau[which(A == j)],
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
    return(list(beta = B[if (burnin != 0){-c(1:burnin)}, ],
                subgrad = S[if(burnin != 0){-c(1:burnin)}, ],
                  acceptHistory = rbind(nAccept, nProp)))
  }
  if (n<p) {
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
    VRCB <- t(VR) %*% C %*% pointEstimate

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
      S0.prop <- lsei(G=G,H=H,E=E,F=F1)
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

      if (!is.null(SFindex) && randomSFindex && (t%%updateSFindex==0) ) {
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
  return(list(beta = B[if (burnin != 0){-c(1:burnin)}, ],
              subgrad = S[if (burnin != 0){-c(1:burnin)}, ],
              pluginbeta = pointEstimate,
              acceptHistory = rbind(nAccept, nProp),
              signchange=nSignChange)) }

}

fixedA.MCMC <- function(X,pointEstimate,sig2,weights=rep(1,max(group)),lbd,group,niter=2000,A,B0,S0,tau,verbose=FALSE) {
  K <- 10
  W <- rep(weights,table(group))
  Psi <- 1/n * crossprod(X)
  if (n > p) {
    inv.Var <- n/sig2 * solve(Psi)
  }

  r.seq <- matrix(, niter, max(group))
  S.seq <- matrix(, niter ,p)
  nAccept <- numeric(2)

  rcur <- group.norm2(B0,group)
  r.seq[1, ] <- rcur
  S.seq[1, ] <- Scur <- S0

  if (n >= p) {
    for (i in 2:niter) {
      r.new <- ld.Update.r(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,tau)
      r.seq[i,] <- rcur <- r.new$r
      Hcur <- r.new$Hcur
      nAccept[1] <- nAccept[1] + r.new$nrUpdate

      S.new <- ld.Update.S (rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,p)
      S.seq[i,] <- Scur <- S.new$S
      Hcur <- S.new$Hcur
      nAccept[2] <- nAccept[2] + S.new$nSUpdate

      if (verbose && (i %% round(niter/10) == 0)) {
        cat("MCMC step,", K, "% Finished\n")
        K <- K+10
      }
    }
  } else {
    for (i in 2:niter) {
      r.new <- hd.Update.r(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,1)
      r.seq[i,] <- rcur <- r.new$r
      Hcur <- r.new$Hcur
      nAccept[1] <- nAccept[1] + r.new$nrUpdate

      S.new <- hd.Update.S (rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,p)
      S.seq[i,] <- Scur <- S.new$S
      Hcur <- S.new$Hcur
      nAccept[2] <- nAccept[2] + S.new$nSUpdate

      if (verbose && (i %% round(niter/10) == 0)) {
        cat("MCMC step,", K, "% Finished\n")
        K <- K+10
      }
    }
  }
  return(list(group.l2.norm = r.seq, subgradient = S.seq, nAccept = nAccept))
}

