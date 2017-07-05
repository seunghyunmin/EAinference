####################################################
# Shrink Y towards zero by Stein estimate
# Y: n*1 matrix
# df: degress of freedom
# return: list(mu: estimate of true mu, L: sure (stein unbiased risk estimate))
ShrinkToZeroStein <- function(Y, df)
{
  B <- df / sum(Y^2)
  mu <- (1 - B) * Y
  L <- max(1 - B, 0)
  return(list(mu = mu, L = L))
}
####################################################
# calculate c(alpha)
# mu: n*1 matrix
# df: degress of freedom
# quantile:
# projection: identity matrix by default or a n*n matrix
# iteNum
# return: the quantile of c(alpha) in Professor's note
quantile_stein <- function(mu, df, quantile, projection=NULL, iteNum = 100)
{
  n <- dim(mu)[1]
  if (is.null(projection))
    projection <- diag(rep(1, n))
  z <- rep(0, iteNum)
  for (i in 1:iteNum)
  {
    Y_prime <- mu + projection %*% matrix(rnorm(n), ncol = 1)
    res <- ShrinkToZeroStein(Y_prime, df)
    z[i] <- res$L - sum((res$mu - mu)^2) / df
  }
  c <- quantile(abs(z) * sqrt(df), probs = quantile)
  return (c)
}
####################################################
# We split X and Y into 2 pieces. (X1, X2) (Y1, Y2).
# Use (X1,Y1) to get the active set A
# Then draw inference using (X2, Y2).

# Given a sequence of ratios, X, Y, lbd, calculate r_s r_w, mu_s mu_w which minimize volumne
# X1, X2: n*p matrx
# Y1, Y2: n*1 matrix
# lbd: tunning parameter, lambda, in lasso
# ratios: a sequence of ratios of tau against lambda (tau the threshold value of \hat{beta})
# alpha:
# iteNum:
# return: list(r_s, r_w, mu_s, mu_w)
projectionStein <- function(X, Y, lbd, ratios, alpha=0.05, iteNum=1000)
{

  W <- sort(sample(1:nrow(X),round(nrow(X)/2)))
  X1 <- X[W, , drop=FALSE]; Y1 <- Y[W, , drop=FALSE]
  X2 <- X[-W, , drop=FALSE]; Y2 <- Y[-W, , drop=FALSE]

  # apply lasso regression of Y onto X
  n <- dim(X2)[1]
  LassoResult <- Lasso.MHLS(X = X1, Y = Y1, lbd = lbd)$B0

  # estimate c by A=0
  stein <- ShrinkToZeroStein(Y2, n)
  c <- quantile_stein(stein$mu, n, 1-alpha/2, iteNum = iteNum)

  # information along ratios
  r_s_ratio <- rep(0, length(ratios))
  r_w_ratio <- rep(0, length(ratios))
  logVol_ratio <- rep(0, length(ratios))

  for (i in 1:length(ratios))
  {
    A <- which(abs(LassoResult) > ratios[i] * lbd)
    if (length(A) != 0)
    {
      r_s_ratio[i] <- sqrt(qchisq(1-alpha/2, df=length(A)) / n)
      P_A <- X2[,A] %*% solve(t(X2[,A]) %*% X2[,A]) %*% t(X2[,A])
      P_perp <- diag(rep(1,n)) - P_A
    } else
    {
      r_s_ratio[i] <- 0
      P_perp <- diag(rep(1,n))
    }
    nMinusA <- n - length(A)

    stein <- ShrinkToZeroStein(P_perp %*% Y2, nMinusA)
    r_w_ratio[i] <- sqrt(nMinusA / n * (c / sqrt(nMinusA) + stein$L))

    #update r_s_ratio r_w_ratio if r_s_ratio exists and calculate volume
    if (length(A) != 0)
    {
      r_s_ratio[i] <- r_s_ratio[i] * sqrt(n / length(A))
      r_w_ratio[i] <- r_w_ratio[i] * sqrt(n / (n - length(A)))
      logVol_ratio[i] <- length(A) * log(r_s_ratio[i]) + nMinusA * log(r_w_ratio[i])
    } else
    {
      logVol_ratio[i] <- nMinusA * log(r_w_ratio[i])
    }
  }
  # information of minimum volume
  i <- which.min(logVol_ratio)
  A <- which(abs(LassoResult) > ratios[i] * lbd)
  # compute mu_s
  if (length(A) != 0)
  {
    P_A <- X2[,A] %*% solve(t(X2[,A]) %*% X2[,A]) %*% t(X2[,A])
    mu_s <- P_A %*% Y2
    P_perp <- diag(rep(1,n)) - P_A
  } else
  {
    mu_s <- 0
    P_perp <- diag(rep(1,n))
  }
  # compute mu_w
  stein <- ShrinkToZeroStein(P_perp %*% Y2, n - length(A))
  mu_w <- stein$mu

  result <- list(r_s=r_s_ratio[i], r_w=r_w_ratio[i], mu_s=mu_s, mu_w=mu_w)
  return(result)
}
