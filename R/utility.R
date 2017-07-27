#' @importFrom gglasso gglasso coef.gglasso
#' @importFrom msm rtnorm
#' @importFrom mvtnorm dmvnorm
#' @import graphics
#' @import stats

# devtools::check()
# devtools::build_win() : to check package under Window enviornment
#-------------------------------------------
# Utility functions for Individual Lasso
#-------------------------------------------
logFlipTR <- function(x,t)
{
  logr <- log(2) + dnorm(x, mean=0, sd=t,log=TRUE);
  # theta(b_j;0,tau_j^2)/(1/2)
  return(logr);
}

SummBeta <- function(x) {
  c( mean=mean(x), median = median(x), s.d = sd(x), quantile(x,c(.025,.975)) )
}

SummSign <- function(x) {
  n=length(x)
  return ( c(Positive.proportion=sum(x>0)/n , Zero.proportion=sum(x==0)/n, Negative.proportion=sum(x<0)/n  ) )
}

UpdateBA <- function(Bcur,Scur,tau,A,I,Rcur,logfRcur,VRC,lbdVRW,InvVarR,
                  tVAN.WA,invB_D,B_F,FlipSA,IndexSF,nA=length(A),nI=length(I))
{
  Bnew=Bcur;
  Bnew[A]=rnorm(nA,Bcur[A],tau)

  diffB=Bnew-Bcur;
  nAcceptBA=0;
  nChangeSign=c(0,0,0);
  for(j in FlipSA)
  {
    if (sign(Bnew[j])!=sign(Bcur[j])) {
      nChangeSign[2]=nChangeSign[2]+1
      Snew=Scur;
      Snew[j]=sign(Bnew[j])
      if (is.null(IndexSF)) {
        Snew[I] = -invB_D%*%tVAN.WA%*%Snew[A]
      } else {
        Snew[I][-IndexSF]=invB_D %*% (-B_F %*% Snew[I][IndexSF] - tVAN.WA%*%Snew[A])
      }

      if (if (is.null(IndexSF)) {all(abs(Snew[I]) <= 1)} else {all(abs(Snew[I][-IndexSF]) <= 1)}) {
        nChangeSign[3]=nChangeSign[3]+1
        diffR=diffB[j]*VRC[,j] + lbdVRW%*%(Snew-Scur);
        Rnew=Rcur+diffR;
        logfRnew=-0.5*sum(Rnew^2*InvVarR);
        logTargetRatio=logfRnew-logfRcur;

        logMH=logTargetRatio;
        u=runif(1);
        if(log(u)<logMH)
        {
          # loc = which(apply(S.temp$Sign,1,function(x) all(x==Scur[A])))
          # if(length(loc)!=0) {
          # 	S.temp$Sign[loc,]=sign(Bcur[A])
          # 	S.temp$Subgrad[loc,]=Scur[I]
          # } else {
          # 	S.temp$Sign=rbind(S.temp$Sign,sign(Bcur[A]))
          # 	S.temp$Subgrad=rbind(S.temp$Subgrad,Scur[I])
          # }
          Rcur=Rnew;
          logfRcur=logfRnew;
          Bcur[j]=Bnew[j];
          Scur=Snew
          nAcceptBA=nAcceptBA+1;
          nChangeSign[1]=nChangeSign[1]+1;
        }
      }
    } else {
      diffR=diffB[j]*VRC[,j];
      Rnew=Rcur+diffR;
      logfRnew=-0.5*sum(Rnew^2*InvVarR);
      logTargetRatio=logfRnew-logfRcur;

      logMH=logTargetRatio;
      u=runif(1);
      if(log(u)<logMH)
      {
        Rcur=Rnew;
        logfRcur=logfRnew;
        Bcur[j]=Bnew[j];
        nAcceptBA=nAcceptBA+1;
      }
    }
  }
  return(list(B=Bcur,S=Scur,Rvec=Rcur,logf=logfRcur,nAccept=nAcceptBA,nChangeSign=nChangeSign));
}

UpdateBA.fixedSA <- function(Bcur,tau,A,Rcur,logfRcur,VRC,InvVarR) # Fix sign(beta[A])
{
  nA=length(A);
  LUbounds=matrix(0,nA,2);
  LUbounds[Bcur[A]>0,2]=Inf;
  LUbounds[Bcur[A]<0,1]=-Inf;
  Bnew=Bcur;
  Bnew[A]=rtnorm(nA,Bcur[A],tau,lower=LUbounds[,1],upper=LUbounds[,2]);

  Ccur=pnorm(0,mean=Bcur[A],sd=tau,lower.tail=TRUE,log.p=FALSE);
  Ccur[Bcur[A]>0]=1-Ccur[Bcur[A]>0];
  Cnew=pnorm(0,mean=Bnew[A],sd=tau,lower.tail=TRUE,log.p=FALSE);
  Cnew[Bcur[A]>0]=1-Cnew[Bcur[A]>0];
  lqratio=log(Ccur/Cnew);

  diffB=Bnew-Bcur;
  i=1;
  nAcceptBA=0;
  for(j in A)
  {
    diffR=diffB[j]*VRC[,j];
    Rnew=Rcur+diffR;
    logfRnew=-0.5*sum(Rnew^2*InvVarR);
    logTargetRatio=logfRnew-logfRcur;

    logMH=logTargetRatio+lqratio[i];
    u=runif(1);
    if(log(u)<logMH)
    {
      Rcur=Rnew;
      logfRcur=logfRnew;
      Bcur[j]=Bnew[j];
      nAcceptBA=nAcceptBA+1;
    }
    i=i+1;
  }
  return(list(B=Bcur,Rvec=Rcur,logf=logfRcur,nAccept=nAcceptBA));
}

# Subfunction for UpdateBA : In high-dim, when sign(BA) changes, update SI.
# ProposeSI=function(PropSA, CurrSI, S.temp, tVAN.WA, nI, E, G, H, iter=10){
# 	loc = which(apply(S.temp$Sign,1,function(x) all(x==PropSA)))
# 	if(length(loc)!=0) return(S.temp$Subgrad[loc,])
# 	# There already is a valid SI
#
# 	F1 <- -tVAN.WA%*%PropSA
#
# 	x0 <- suppressWarnings(lsei(A=diag(nI),B=c(CurrSI),G=G,H=H,E=E,F=F1))
# 	if (x0$IsError == TRUE) return(rep(9,nI)) # Solution does not exist
# 	#return(x0$X[1:nI])
# 	#return(suppressWarnings(xsample(A=diag(nI),B=c(CurrSI),G=G,H=H,E=E,F=F1,iter=iter,x0=x0$X))$X[100,1:nI])
# 	return(suppressWarnings(xsample(G=G,H=H,E=E,F=F1,iter=iter,x0=x0$X))$X[iter,1:nI])
# 	# Draw valid SI
# }

UpdateSI <- function(Scur,A,I,Rcur,n,p,logfRcur,lbdVRW,InvVarR,
                  tVAN.WA,invB_D,BDF, B_F, IndexSF, SIscale=1, nA=length(A))
{
  if (SIscale <1) stop("SIscale has to be greater than 1")

  nAcceptSI=0;
  U=-tVAN.WA%*%Scur[A]
  #V_AN=V[A,N,drop=FALSE];
  #V_IN=V[I,N,drop=FALSE];
  #V_AR=V[A,R,drop=FALSE];
  #V_IR=V[I,R,drop=FALSE];
  #t(V_AN)%*%W[A,A]%*%Snew[A]+t(V_IN)%*%W[I,I]%*%Snew[I]

  a=invB_D%*%U + rep(-1,p-n)
  b=a + rep(2,p-n)

  for ( i in 1:(n-nA)) {

    Scur_F=Scur[I][IndexSF]
    Scur_D=setdiff(Scur[I],Scur_F)

    bd=cbind(a-BDF[,-i,drop=FALSE]%*%Scur_F[-i], b-BDF[,-i,drop=FALSE]%*%Scur_F[-i])
    #		BDF[,i]*Scur_F[i]
    bd=bd/BDF[,i]

    bd=t(apply(bd,1,sort))

    lowbd=max(c(bd[,1],-1))
    highbd=min(c(bd[,2],1))

    if(lowbd>highbd) {
      if ( lowbd-highbd > 1e-10) {
        stop("Lowerbound is greater than Uppderbound?")
      } else {
        lowbd = highbd
      }
    }

    lnth=(highbd-lowbd)/SIscale

    Snew_F=Scur_F
    Snew_F[i]=runif(1,(lowbd+highbd)/2-lnth/2,(lowbd+highbd)/2+lnth/2)
    Snew_D = invB_D%*%(U-B_F%*%Snew_F)

    Snew_I=Scur[I]
    Snew_I[IndexSF]=Snew_F
    Snew_I[-IndexSF]=Snew_D

    diffSI=Snew_I-Scur[I]

    diffR=lbdVRW[,I,drop=FALSE]%*%diffSI;
    Rnew=Rcur+diffR;
    logfRnew=-0.5*sum(Rnew^2*InvVarR);
    logTargetRatio=logfRnew-logfRcur;

    #logPropRatio = log(lnth) -log(2)
    # Proposal density should be the same in both directions.

    #logMH=logTargetRatio + logPropRatio;
    logMH=logTargetRatio;

    u=runif(1);
    if(log(u)<logMH)
    {
      Rcur=Rnew;
      logfRcur=logfRnew;
      Scur[I]=Snew_I;
      nAcceptSI=nAcceptSI+1;
    }
  }
  return(list(S=Scur,Rvec=Rcur,logf=logfRcur,nAccept=nAcceptSI));
}


CalTmat <- function(p,n,V,LBD,W,lbd,R,N,A,I)
{
  V_IN=V[I,N,drop=FALSE];
  V_AR=V[A,R,drop=FALSE];
  V_IR=V[I,R,drop=FALSE];
  if(length(A)<n)
  {
    BNul=as.matrix(svd(t(V_IN)%*%as.matrix(W[I,I]),nv=length(I))$v[,(p-n+1):length(I)]);
    Tmat=cbind(LBD%*%t(V_AR),lbd*t(V_IR)%*%as.matrix(W[I,I])%*%BNul);
    #LBD%*%t(V_AR) == t(VR)%*%C[,A]
  }else{
    BNul=NULL;
    Tmat=LBD%*%t(V_AR);
  }
  logdetT=determinant(as.matrix(Tmat),logarithm=TRUE)$modulus[1];
  return(list(Basisset=BNul,T=Tmat,logdetT=logdetT));
}
#-------------------------------------------
# Utility functions for Group Lasso
#-------------------------------------------
group.norm <- function(x, group, al = 2) {
  tapply(x,group, function(x) sum(abs(x)^al)^(1/al))
}

group.norm2 <- function(x, group) {
  result <- c()
  for (i in 1:max(group)) {
    result[i] <- sqrt(crossprod(x[group==i]))
  }
  return(result)
}

# T(s,A) : p x (n-|A|) matrix s.t. ds = T(s,A)ds_F
TsA <- function(Q, s, group, A, n, p) {
  # even if length(A) == 0, everything will work just fine !!
  # when lengthI(A) == n, we only compute F2 function.
  if (n < p && missing(Q)) {
    stop("High dimensional setting needs Q")
  }
  nA <- length(A)
  rankX <- min(c(n,p))
  Subgradient.group.matix <- matrix(0, length(unique(group)), p)

  for (i in 1:length(unique(group))) {
    Subgradient.group.matix[i, which(group == i)] = s[group == i]
  }
  Subgradient.group.matix <- Subgradient.group.matix[A, ,drop=FALSE]

  if (n >= p) {
    B <- Subgradient.group.matix
  } else {
    B <- rbind(t(Q), Subgradient.group.matix)
  }

  if (nA != 0 ) {
    P <- matrix(0, p, p)
    Permute.order <- 1:p

    if (n < p) {nrowB <- p - n + nA} else
    {nrowB <- nA}
    for (i in 1:nrowB) {
      if (B[i, rankX - nA + i] == 0) {
        W1 <- which(B[i,] !=0)[1]
        #Permute.order[c(W1,n-length(A)+i)] = Permute.order[c(n-length(A)+i,W1)]
        Permute.order[c(which(Permute.order == W1),rankX-nA+i)] =
          Permute.order[c(rankX-nA+i,which(Permute.order == W1))]
        #print(Permute.order)
      }
    }
    for ( i in 1:p) {
      P[Permute.order[i], i] <- 1
    }
    B <- (B%*%P)
  } else {
    P <- diag(p)
  }

  if (nA == 0 && n >= p) return(P);

  if (rankX-nA >= 1) {
    BF <- B[, 1:(rankX - nA),drop=FALSE]
    BD <- B[, -c(1:(rankX - nA)),drop=FALSE]
    #Result <- P %*% rbind(diag(n-length(A)), -solve(BD)%*%BF)
    Result <- P %*% rbind(diag(rankX-nA), -solve(BD,BF))
  } else {
    Result <- -solve(B)
  }
  return(Result)
}
#------------------------------
# TsA.qr is not updated
# RD is not invertible-guaranteed.
TsA.qr <- function(Q, s, group, A, n, p) { # T(s,A)
  # even if length(A) == 0, everything will work just fine !!
  # when length(A) == n, we only compute F2 function.
  if (n == length(A)) {stop("|A| should be smaller than n")}
  if (n < p && missing(Q)) {
    stop("High dimensional setting needs Q")
  }

  nA <- length(A)

  Subgradient.group.matix <- matrix(0, length(unique(group)), p)
  for (i in 1:length(unique(group))) {
    Subgradient.group.matix[i, which(group == i)] = s[group == i]
  }
  Subgradient.group.matix <- Subgradient.group.matix[A, , drop=FALSE]
  if (n >= p) {
    B <- Subgradient.group.matix
  } else {
    B <- rbind(t(Q), Subgradient.group.matix)
  }

  if (n < p) { nrowB <- p - n + nA ; rankX <- n} else
  {nrowB <- nA; rankX <- p}

  QR.B <- qr(B)
  Pivot <- sort.list(qr(B)$pivot)
  B.Q <- qr.Q(QR.B)
  #B.R <- qr.R(QR.B)[,Pivot]
  B.R <- qr.R(QR.B)
  tP <- matrix(0, p, p)

  for (i in 1:p) {
    tP[Pivot[i], i] = 1
  }

  P <- t(tP)

  if (rankX-nA >= 1) {
    RF <- B.R[, 1:(rankX - nA), drop=FALSE]
    RD <- B.R[, -c(1:(rankX - nA)), drop=FALSE]
    Result <- P %*% rbind(diag(rankX-nA), -solve(RD)%*%RF)
  }
  return(Result)
}

# TsA.null <- function(t.XWinv, s, group, A, n, p) { # T(s,A)
#   # even if length(A) == 0, everything will work just fine !!
#   # when length(A) == n, we only compute F2 function.
#   # Updated for Low-dim case
#   if (missing(t.XWinv) && n < p) {
#     stop("When n < p, t.XWinv is needed")
#   }
#
#   nA <- length(A)
#   if (n < p) { # High-dim
#     if (nA !=0) {
#       Subgradient.group.matix <- matrix(0, nA, p)
#       for (i in 1:nA) {
#         Subgradient.group.matix[i, which(group == A[i])] <- s[group == A[i]]
#       }
#       t.XWinv %*% Null(t(Subgradient.group.matix %*% t.XWinv))
#     } else {
#       t.XWinv
#     }
#   } else { # Low-dim
#     if (nA !=0) {
#       Subgradient.group.matix <- matrix(0, nA, p)
#       for (i in 1:nA) {
#         Subgradient.group.matix[i, which(group == A[i])] <- s[group == A[i]]
#       }
#       Null(t(Subgradient.group.matix))
#     } else {
#       diag(p)   # if |A| = 0, T should be pXp identity matrix.
#     }
#   }
# }

# F1 = r \circ \psi , eq(3.6), p x p matrix
F1 <- function(r, Psi, group) {
  Result <- Psi
  for (i in 1:length(unique(group))) {
    Result[, group == i] <- Result[, group == i] * r[i]
  }
  return(Result)
}

# F2 = \psi \circ \eta , eq(3.7), p x J matrix, where J is the number of groups
F2 <- function(s, Psi, group) {
  Result <- matrix(, nrow(Psi), length(unique(group)))
  for (i in 1:length(unique(group))) {
    Result[, i] = crossprod(t(Psi[, group == i]), s[group == i])
  }
  return(Result)
}

log.Jacobi.partial <- function(X, s, r, Psi, group, A, lam, W, TSA) { # log(abs(det(X %*% [F2[,A] | (F1 + lam * W) %*% TsA])))
  n <- nrow(X)
  p <- ncol(X)
  table.group <- table(group)

  # W <- c()
  # for (i in 1:length(table.group)) {
  #   W <- c(W, rep(weights[i], table.group[i]))
  # }

  if (nrow(X) < ncol(X)) { # High-dim
    if (n == length(A)) {
      log.Det <- determinant(X %*% F2(s, Psi, group)[,A])
    } else {
      log.Det <- determinant(X %*% cbind(F2(s, Psi, group)[, A], (F1(r, Psi, group) + lam * diag(W)) %*% TSA))
    }
    return(log.Det[[1]][1]);
  } else { # Low-dim
    if (p == length(A)) {
      log.Det <- determinant(F2(s,Psi,group))
    } else {
      log.Det <- determinant(cbind(F2(s, Psi, group)[, A], (F1(r, Psi, group) + lam * diag(W)) %*% TSA))
    }
    return(log.Det[[1]][1]);
  }
}

ld.Update.r <- function(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,tau,PEtype,n,p) {
  rprop <- rcur;
  nrUpdate <- 0;
  Bcur <- Bprop <- Scur * rep(rcur,table(group));
  TSA.cur <- TSA.prop <- TsA(,Scur,group,A,n,p);
  for (i in A) {
      #rprop[i] <- truncnorm::rtruncnorm(1, 0, , rcur[i], sqrt(tau[which(A==i)] * ifelse(rcur[i]!=0,rcur[i],1)))
    rprop[i] <- rtnorm(n = 1, mean = rcur[i], sd = sqrt(tau[which(A==i)] * ifelse(rcur[i]!=0,rcur[i],1)), lower = 0, upper = Inf)
    Bprop[group==i] <- rprop[i] * Scur[group==i]

    if (PEtype == "coeff") {
      Hprop <- drop(Psi %*% drop(Bprop - pointEstimate) + lbd * W * drop(Scur))
    } else {
      Hprop <- drop(Psi %*% drop(Bprop) - t(X) %*% pointEstimate / n + lbd * W * drop(Scur))
    }

    Hdiff <- Hcur - Hprop

    lNormalRatio <- drop(t(Hdiff)%*% inv.Var %*% (Hprop + Hdiff/2))
    #dmvnorm(Hprop,,sig2 / n * Psi,log=T) - dmvnorm(Hcur,,sig2 / n * Psi,log=T)
    lJacobianRatio <- log.Jacobi.partial(X,Scur,rprop,Psi,group,A,lbd,W,TSA.prop) -
      log.Jacobi.partial(X,Scur,rcur,Psi,group,A,lbd,W,TSA.cur)
    lProposalRatio <- pnorm(0,rcur[i],sqrt(tau[which(A==i)] * rcur[i]), lower.tail=FALSE, log.p=TRUE) -
      pnorm(0,rprop[i],sqrt(tau[which(A==i)] * rprop[i]), lower.tail=FALSE, log.p=TRUE)
    lAcceptanceRatio <-  lNormalRatio + lJacobianRatio + lProposalRatio
    if (lAcceptanceRatio <= log(runif(1))) { # Reject
      rprop[i] <- rcur[i];
      Bprop[group==i] <- Bcur[group==i]
    } else { # Accept
      nrUpdate <- nrUpdate + 1;
      Hcur <- Hprop;
    }
  }
  return(list(r = rprop, Hcur = Hcur, nrUpdate = nrUpdate))
}
ld.Update.S <- function(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,PEtype,n,p) {
  Sprop <- Scur;
  nSUpdate <- 0;
  #p <- ncol(X)
  for (i in 1:max(group)) {
    if (i %in% A) {Sprop[group == i] <- rUnitBall.surface(sum(group == i))} else {
      Sprop[group ==i] <- rUnitBall(sum(group==i))
    }

    if (PEtype == "coeff") {
      Hprop <- drop(Psi %*% drop(Sprop * rep(rcur,table(group)) - pointEstimate) + lbd * W * drop(Sprop))
    } else {
      Hprop <- drop(Psi %*% drop(Sprop * rep(rcur,table(group))) - t(X) %*% pointEstimate / n + lbd * W * drop(Sprop))
    }
    Hdiff <- Hcur - Hprop

    lNormalRatio <- drop(t(Hdiff)%*% inv.Var %*% (Hprop + Hdiff/2))
    #dmvnorm(Hprop,,sig2 / n * Psi,log=T) - dmvnorm(Hcur,,sig2 / n * Psi,log=T)
    lJacobianRatio <- log.Jacobi.partial(X,Sprop,rcur,Psi,group,A,lbd,W,TsA(,Sprop,group,A,n,p)) -
      log.Jacobi.partial(X,Scur,rcur,Psi,group,A,lbd,W,TsA(,Scur,group,A,n,p))
    lAcceptanceRatio <-  lNormalRatio + lJacobianRatio
    if (lAcceptanceRatio <= log(runif(1))) { # Reject
      Sprop[group == i] <- Scur[group == i];
    } else { # Accept
      nSUpdate <- nSUpdate + 1;
      Hcur <- Hprop;
    }
  }
  return(list(S = Sprop, Hcur = Hcur, nSUpdate = nSUpdate))
}
rUnitBall.surface <- function(p) {
  x <- rnorm(p)
  x / sqrt(crossprod(x))
}
rUnitBall <- function(p) {
  x <- rnorm(p,,1/sqrt(2));
  y <- rexp(1)
  x / sqrt(y+crossprod(x))
}
# hd.Update.r <- function(rcur,Scur,A,Hcur,X,coeff,Psi,W,lbd,group,inv.Var,tau) {}
# hd.Update.S <- function(rcur,Scur,A,Hcur,X,coeff,Psi,W,lbd,group,inv.Var,p) {}
#-------------------------------------------
# Utility functions for MHLS summary
#-------------------------------------------
SummBeta <- function ( x ) {
  c( mean=mean(x) , median = median(x) , s.d = sd(x) , quantile(x,c(.025,.975)) )
}

SummSign <- function ( x ) {
  n=length(x)
  return ( c(Positive.proportion=sum(x>0)/n , Zero.proportion=sum(x==0)/n, Negative.proportion=sum(x<0)/n  ) )
}
#-------------------------------------------
# Utility functions for scaled lasso / scaled group lasso
#-------------------------------------------
TsA.slasso <- function(SVD.temp, Q, s, W, group, A, n, p) {
  # even if length(A) == 0, everything will work just fine !!
  # when lengthI(A) == n, we only compute F2 function.
  if (n < p && missing(Q)) {
    stop("High dimensional setting needs Q")
  }
  nA <- length(A)
  Subgradient.group.matix <- matrix(0, nA, p)

  IndWeights <- W
  if (nA != 0) {
    for (i in 1:nA) {
      Subgradient.group.matix[i, which(group == A[i])] = s[group == A[i]]
    }
  }
  #Subgradient.group.matix <- Subgradient.group.matix[A, ,drop=FALSE]

  #all.equal(V%*%diag(1/D^2)%*%t(V) , t(X) %*% solve(tcrossprod(X)) %*% solve(tcrossprod(X)) %*% X)
  #all.equal(U%*%diag(D)%*%t(V) , X)

  B <- rbind(t(Q), Subgradient.group.matix,
             t(s * IndWeights) %*% SVD.temp %*% diag(IndWeights))

  if (nA != 0 ) {
    P <- matrix(0, p, p)
    Permute.order <- 1:p

    # if (n < p) {nrowB <- p - n + nA} else
    # {nrowB <- nA}
    for (i in 1:nrow(B)) {
      if (B[i, n - nA - 1 + i] == 0) {
        W1 <- which(B[i,] !=0)[1]
        #Permute.order[c(W1,n-length(A)+i)] = Permute.order[c(n-length(A)+i,W1)]
        Permute.order[c(which(Permute.order == W1),n-nA-1+i)] =
          Permute.order[c(n-nA-1+i,which(Permute.order == W1))]
        #print(Permute.order)
      }
    }
    for ( i in 1:p) {
      P[Permute.order[i], i] <- 1
    }
    B <- (B%*%P)
  } else {
    P <- diag(p)
  }

  if (n-nA-1 >= 1) { # n-nA-1 : # of free coordinate
    BF <- B[, 1:(n - nA - 1),drop=FALSE]
    BD <- B[, -c(1:(n - nA - 1)),drop=FALSE]
    #Result <- P %*% rbind(diag(n-length(A)), -solve(BD)%*%BF)
    Result <- P %*% rbind(diag(n-nA-1), -solve(BD,BF))
  } else {
    Result <- -solve(B)
  }
  return(Result)
}
log.Jacobi.partial.slasso <- function(X, s, r, Psi, group, A, lam, hsigma, W, TSA) { # log(abs(det(X %*% [F2[,A] | (F1 + lam * W) %*% TsA])))
  # This function is only for high-dimensional cases.
  n <- nrow(X)
  p <- ncol(X)
  table.group <- table(group)
  if (n > p) stop("High dimensional setting is required.")
  if (n == (length(A)-1)) {
    log.Det <- determinant(X %*% cbind(F2(s, Psi, group)[,A], diag(lam * hsigma * W) %*% s))
  } else {
    log.Det <- determinant(X %*% cbind(F2(s, Psi, group)[, A], (F1(r, Psi, group)
                                                                + diag(lam * hsigma * W)) %*% TSA,diag(lam * W) %*% s))
  }
  return(log.Det[[1]][1]);
}

#-------------------------------------------
# Utility functions for hdi package.
#-------------------------------------------

## This file contains functions that calculate the Z vectors as nodewise lasso
## residuals
## see http://arxiv.org/abs/1110.2563
## and can also calculate the Thetahat matrix as in
## http://arxiv.org/abs/1303.0518

score.nodewiselasso <- function(x,
                                wantTheta = FALSE,
                                verbose = FALSE,
                                lambdaseq = "quantile",
                                parallel = FALSE,
                                ncores = 8,
                                oldschool = FALSE,
                                lambdatuningfactor = 1,
                                cv.verbose = FALSE)
{
  ## Purpose:
  ## This function calculates the score vectors Z OR the matrix of nodewise
  ## regression coefficients Thetahat, depending on the argument wantTheta.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  ## First, a sequence of lambda values over all nodewise regressions is created

  lambdas <-
    switch(lambdaseq,
           "quantile" = nodewise.getlambdasequence(x), ## this is preferred
           "linear"   = nodewise.getlambdasequence.old(x,verbose),
           ## otherwise
           stop("invalid 'lambdaseq': ", lambdaseq))

  if(verbose){
    cat("Using the following lambda values:", lambdas, "\n")
  }

  ## 10-fold cv is done over all nodewise regressions to calculate the error
  ## for the different lambda
  cvlambdas <- cv.nodewise.bestlambda(lambdas = lambdas, x = x,
                                      parallel = parallel, ncores = ncores,
                                      oldschool = oldschool,
                                      verbose = cv.verbose)
  if(verbose){
    cat(paste("lambda.min is", cvlambdas$lambda.min), "\n")
    cat(paste("lambda.1se is", cvlambdas$lambda.1se), "\n")
  }

  if(lambdatuningfactor == "lambda.1se"){
    if(verbose)
      cat("lambda.1se used for nodewise tuning\n")
    ## We use lambda.1se for bestlambda now!!!
    bestlambda <- cvlambdas$lambda.1se
  }else{
    if(verbose)
      cat("lambdatuningfactor used is", lambdatuningfactor, "\n")

    bestlambda <- cvlambdas$lambda.min * lambdatuningfactor
  }

  if(verbose){
    cat("Picked the best lambda:", bestlambda, "\n")
    ##print("with the error ")
    ##print(min(err))
  }

  ## Having picked the 'best lambda', we now generate the final Z or Thetahat
  if(wantTheta){
    out <- score.getThetaforlambda(x = x,
                                   lambda = bestlambda,
                                   parallel = parallel,
                                   ncores = ncores,
                                   oldschool = TRUE,
                                   verbose = verbose)
  }else{
    Z <- score.getZforlambda(x = x, lambda = bestlambda, parallel = parallel,
                             ncores = ncores, oldschool = oldschool)
    out <- Z
  }
  return.out <- list(out = out,
                     bestlambda = bestlambda)
  return(return.out)
}

score.getThetaforlambda <- function(x, lambda, parallel = FALSE, ncores = 8,
                                    oldschool = FALSE, verbose = FALSE,
                                    oldtausq = TRUE)
{
  ## Purpose:
  ## This function is for calculating Thetahat once the desired tuning
  ## parameter is known
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  message("Calculating Thetahat by doing nodewise regressions and dropping the unpenalized intercept")
  n <- nrow(x)
  p <- ncol(x)
  C <- diag(rep(1,p))
  T2 <- numeric(p)

  if(oldschool){
    message("doing getThetaforlambda oldschool")
    for(i in 1:p){
      glmnetfit <- glmnet(x[,-i], x[,i])
      coeffs <- as.vector(predict(glmnetfit,x[,-i], type = "coefficients",
                                  s = lambda))[-1]
      ## we just leave out the intercept

      C[-i,i] <- -as.vector(coeffs)
      if(oldtausq){
        ## possibly quite incorrect,it ignores the intercept, but intercept is
        ## small anyways. Initial simulations show little difference.
        T2[i] <- as.numeric(crossprod(x[,i]) / n - x[,i] %*% (x[,-i] %*%
                                                                coeffs) / n)
      }else{
        ##print("now doing the new way of calculating tau^2")
        T2[i] <- as.numeric((x[,i] %*%
                               (x[,i] - predict(glmnetfit,x[,-i],s =lambda)))/n)
      }
    }
  }else{
    stop("not implemented yet!")
  }
  thetahat <- C %*% solve(diag(T2))
  if(verbose){
    cat("1/tau_j^2:", solve(diag(T2)), "\n")
  }
  ##this is thetahat ^ T!!
  thetahat <- t(thetahat)

  if(all(thetahat[lower.tri(thetahat)] == 0,
         thetahat[upper.tri(thetahat)] == 0) && verbose)
    cat("Thetahat is a diagonal matrix!\n")

  return(thetahat)
}

score.getZforlambda <- function(x, lambda, parallel = FALSE, ncores = 8,
                                oldschool = FALSE)
{
  ## Purpose:
  ## This function is for calculating Z once the desired tuning parameter is
  ## known
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version)

  n <- nrow(x)
  p <- ncol(x)
  Z <- matrix(numeric(n*p),n)

  if(oldschool){
    message("doing getZforlambda oldschool")
    for(i in 1:p){
      glmnetfit <- glmnet(x[,-i],x[,i])
      prediction <- predict(glmnetfit,x[,-i],s=lambda)
      Z[,i] <- x[,i] - prediction
    }
  }else{
    ## REPLACING THE FOR LOOP
    if(parallel){
      Z <- mcmapply(score.getZforlambda.unitfunction, i = 1:p, x = list(x = x),
                    lambda = lambda, mc.cores = ncores)

    }else{
      Z <- mapply(score.getZforlambda.unitfunction, i = 1:p, x = list(x = x),
                  lambda = lambda)
    }
  }
  ## rescale Z such that t(Zj) Xj/n = 1 \-/ j
  Z <- score.rescale(Z,x)
  return(Z)
}

score.getZforlambda.unitfunction <- function(i, x, lambda)
{
  ## Purpose:
  ## Calculate the residuals of a nodewise regression of one column vs the
  ## other columns
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  glmnetfit  <- glmnet(x[,-i],x[,i])
  prediction <- predict(glmnetfit,x[,-i],s=lambda)
  return(x[,i] - prediction)
}

score.rescale <- function(Z, x)
{
  ## Purpose:
  ## Rescale the Z appropriately such that such that t(Zj) Xj/n = 1 for all j
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  scaleZ <- diag(crossprod(Z,x)) / nrow(x)
  Z      <- scale(Z, center = FALSE, scale = scaleZ)
  return(list(Z=Z,scaleZ=scaleZ))
}

nodewise.getlambdasequence <- function(x)
{
  ## Purpose:
  ## this method returns a lambdasequence for the nodewise regressions
  ## by looking at the automatically selected lambda sequences
  ## for each nodewise regression by glmnet.
  ## Equidistant quantiles of the complete set of lambda values are returned.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  nlambda <- 100 ## use the same value as the glmnet default
  p <- ncol(x)

  lambdas <- c()
  for(c in 1:p){
    lambdas <- c(lambdas,glmnet(x[,-c],x[,c])$lambda)
  }

  lambdas <- quantile(lambdas, probs = seq(0, 1, length.out = nlambda))
  lambdas <- sort(lambdas, decreasing = TRUE)
  return(lambdas)
}

cv.nodewise.err.unitfunction <- function(c, K, dataselects, x, lambdas,
                                         verbose, p) {
  ## Purpose:
  ## this method returns the K-fold cv error made by the nodewise regression
  ## of the single column c of x vs the other columns for all values of the
  ## tuning parameters in lambdas.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  if(verbose){ ##print some information out about the progress
    ##report every 25%
    interesting.points <- round(c(1/4,2/4,3/4,4/4)*p)
    names(interesting.points) <- c("25%","50%","75%","100%")
    if(c %in% interesting.points){
      cat("The expensive computation is now",
          names(interesting.points)[c == interesting.points],
          "done\n")
    }
  }

  ## return  'totalerr'
  cv.nodewise.totalerr(c = c,
                       K = K,
                       dataselects = dataselects,
                       x = x,
                       lambdas = lambdas)
}

## gets the standard error for one particular lambda
cv.nodewise.stderr <- function(K, x, dataselects, lambda, parallel, ncores)
{
  ## Purpose:
  ## this method returns the standard error
  ## of the average K-fold cv error made by the nodewise regression
  ## of each column vs the other columns for a single lambda value.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  p <- ncol(x)
  if(parallel){
    totalerr <- mcmapply(cv.nodewise.totalerr,
                         c = 1:p,
                         K = K,
                         dataselects = list(dataselects = dataselects),
                         x = list(x = x),
                         lambdas = lambda,
                         mc.cores = ncores)
  }else{
    totalerr <- mapply(cv.nodewise.totalerr,
                       c = 1:p,
                       K = K,
                       dataselects = list(dataselects = dataselects),
                       x = list(x = x),
                       lambdas = lambda)
  }
  ## get the mean over the variables
  totalerr.varmean <- rowMeans(totalerr)

  ## get the stderror over the K;  return stderr.forlambda
  sd(totalerr.varmean) / sqrt(K)
}

cv.nodewise.totalerr <- function(c, K, dataselects, x, lambdas)
{
  ## Purpose:
  ## this method returns the error made for each fold of a K-fold cv
  ## of the nodewise regression of the single column c of x vs the other
  ## columns for all values of the tuning parameters in lambdas.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  totalerr <- matrix(nrow = length(lambdas), ncol = K)

  for(i in 1:K){ ## loop over the test sets
    whichj <- dataselects == i ##the test part of the data

    glmnetfit <- glmnet(x = x[!whichj,-c, drop = FALSE],
                        y = x[!whichj, c, drop = FALSE],
                        lambda = lambdas)
    predictions  <- predict(glmnetfit, newx = x[whichj, -c, drop = FALSE],
                            s = lambdas)
    totalerr[, i] <- apply((x[whichj, c] - predictions)^2, 2, mean)
  }

  totalerr
}


cv.nodewise.bestlambda <- function(lambdas, x, K = 10, parallel = FALSE,
                                   ncores = 8, oldschool = FALSE,
                                   verbose = FALSE)
{
  ## Purpose:
  ## this function finds the optimal tuning parameter value for minimizing
  ## the K-fold cv error of the nodewise regressions.
  ## A second value of the tuning parameter, always bigger or equal to the
  ## former, is returned which is calculated by allowing the cv error to
  ## increase by the amount of
  ## 1 standard error (a similar concept as to what is done in cv.glmnet).
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),

  n <- nrow(x)
  p <- ncol(x)
  l <- length(lambdas)

  ## Based on code from cv.glmnet for sampling the data
  dataselects <- sample(rep(1:K, length = n))

  if(oldschool){
    message("doing cv.nodewise.error oldschool")
    totalerr <- numeric(l)
    for(c in 1:p){ ## loop over the nodewise regressions
      for(i in 1:K){ ## loop over the test sets
        whichj <- dataselects == i ## the test part of the data

        glmnetfit <- glmnet(x[!whichj,-c,drop=FALSE], x[!whichj,c,drop=FALSE],
                            lambda=lambdas)
        predictions <- predict(glmnetfit,x[whichj, -c, drop = FALSE],
                               s = lambdas)
        totalerr <- totalerr + apply((x[whichj,c]-predictions)^2, 2, mean)
      }
    }
    totalerr <- totalerr / (K * p)
    stop("lambda.1se not implemented for oldschool cv.nodewise.bestlamba")
  }else{
    ## REPLACING THE FOR LOOP

    ##totalerr <- matrix(nrow = l, ncol = p)

    if(parallel){
      totalerr <- mcmapply(cv.nodewise.err.unitfunction,
                           c = 1:p,
                           K = K,
                           dataselects = list(dataselects = dataselects),
                           x = list(x = x),
                           lambdas = list(lambdas = lambdas),
                           mc.cores = ncores,
                           SIMPLIFY = FALSE,
                           verbose = verbose,
                           p=p)
    }else{
      totalerr <- mapply(cv.nodewise.err.unitfunction,
                         c = 1:p,
                         K = K,
                         dataselects = list(dataselects = dataselects),
                         x = list(x = x),
                         lambdas = list(lambdas = lambdas),
                         SIMPLIFY = FALSE,
                         verbose = verbose,
                         p = p)
    }
    ## Convert into suitable array (lambda, cv-fold, predictor)
    err.array  <- array(unlist(totalerr), dim = c(length(lambdas), K, p))
    err.mean   <- apply(err.array, 1, mean) ## 1 mean for each lambda

    ## calculate mean for every lambda x fold combination (= average over p)
    ## for every lambda then get the standard errors (over folds)
    err.se     <- apply(apply(err.array, c(1, 2), mean), 1, sd) / sqrt(K)
    ##totalerr <- apply(totalerr, 1, mean)
  }

  pos.min    <- which.min(err.mean)
  lambda.min <- lambdas[pos.min]

  stderr.lambda.min <- err.se[pos.min]
  ##-   stderr.lambda.min <- cv.nodewise.stderr(K = K,
  ##-                                           x = x,
  ##-                                           dataselects = dataselects,
  ##-                                           lambda = lambda.min,
  ##-                                           parallel = parallel,
  ##-                                           ncores = ncores)
  list(lambda.min = lambda.min,
       lambda.1se = max(lambdas[err.mean < (min(err.mean) + stderr.lambda.min)]))
}

##DEPRECATED, not really the best choice
nodewise.getlambdasequence.old <- function(x,verbose=FALSE)
{
  ## Purpose:
  ## this method returns a lambdasequence for the nodewise regressions
  ## by looking at the automatically selected lambda sequences
  ## for each nodewise regression by glmnet.
  ## It returns a __linear__ interpolation of lambda values between the max and
  ## min lambda value found.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 27 Nov 2012 (initial version),
  nlambda <- 100#take the same value as glmnet does automatically
  p <- ncol(x)
  maxlambda <- 0
  minlambda <- 100

  for(c in 1:p){
    lambdas <- glmnet(x[,-c],x[,c])$lambda

    ##DEBUG
    if(verbose || is.nan(max(lambdas))){
      cat(paste("c:", c, "\n"))
      cat("lambdas:", lambdas, "\n")
      cat("max(lambdas) max(lambdas,na.rm=TRUE) maxlambda: ",
          max(lambdas), max(lambdas,na.rm=TRUE), maxlambda, "\n")
    }
    if(max(lambdas,na.rm=TRUE) > maxlambda){
      maxlambda <- max(lambdas,na.rm=TRUE)
    }
    if(min(lambdas,na.rm=TRUE) < minlambda){
      minlambda <- min(lambdas, na.rm = TRUE)
    }
  }

  lambdas <- seq(minlambda, maxlambda, by = (maxlambda-minlambda)/nlambda)
  ## return
  sort(lambdas, decreasing=TRUE)
}
