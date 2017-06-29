#' @importFrom gglasso gglasso coef.gglasso

# devtools:load_all()
# Cmd + Shift + L , to load the package.
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
    # QQ : LBD%*%t(V_AR)?
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

TsA <- function(Q, s, group, A, n, p) { # T(s,A)
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

F1 <- function(r, Psi, group) { # r \circ \psi , eq(3.6)
  Result <- Psi
  for (i in 1:length(unique(group))) {
    Result[, group == i] <- Result[, group == i] * r[i]
  }
  return(Result)
}

F2 <- function(s, Psi, group) { # \psi \circ \eta , eq(3.7)
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

ld.Update.r <- function(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,tau,type) {
  rprop <- rcur;
  nrUpdate <- 0;
  Bcur <- Bprop <- Scur * rep(rcur,table(group));
  TSA.cur <- TSA.prop <- TsA(,Scur,group,A,n,p);
  for (i in A) {
    rprop[i] <- truncnorm::rtruncnorm(1, 0, , rcur[i], sqrt(tau[which(A==i)] * ifelse(rcur[i]!=0,rcur[i],1)))
    Bprop[group==i] <- rprop[i] * Scur[group==i]

    if (type == "coeff") {
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

ld.Update.S <- function(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,p,type) {
  Sprop <- Scur;
  nSUpdate <- 0;
  #p <- ncol(X)
  for (i in 1:max(group)) {
    if (i %in% A) {Sprop[group == i] <- rUnitBall.surface(sum(group == i))} else {
      Sprop[group ==i] <- rUnitBall(sum(group==i))
    }

    if (type == "coeff") {
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
#=========================================
Test.stats <- function(beta, group) {
  Group.norm <- group.norm2(beta,group)
  Result <- c(sum(Group.norm), max(Group.norm), Group.norm[1], Group.norm[2])
  #names(Result) <- paste("T", 1:4, sep = "")
  names(Result) <- c("sum.group.norm","max.group.norm","1st.group.norm","2nd.group.norm")
  return(Result)
}
#=========================================
Test.stats.A <- function(beta, group, A) {
  Group.norm <- group.norm2(beta,group)
  Result <- c(sum(Group.norm), max(Group.norm), Group.norm[A])
  #names(Result) <- paste("T", 1:4, sep = "")
  names(Result) <- c("sum.group.norm","max.group.norm",paste("group.norm",A,sep=""))
  return(Result)
}

#====================================
# MHLS class
#====================================
SummBeta <- function ( x ) {
  c( mean=mean(x) , median = median(x) , s.d = sd(x) , quantile(x,c(.025,.975)) )
}

SummSign <- function ( x ) {
  n=length(x)
  return ( c(Positive.proportion=sum(x>0)/n , Zero.proportion=sum(x==0)/n, Negative.proportion=sum(x<0)/n  ) )
}
