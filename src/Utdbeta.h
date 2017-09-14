#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::rowvec Utdbeta(NumericMatrix X, NumericVector Y, NumericMatrix XY, arma::rowvec Beta)
{
  int n = Y.size();
  int p = X.ncol();
  NumericVector TEMP(p);

  for (int i=0; i<n; i++) {
    double xb = 0;
    for (int j=0; j<p; j++) {
      xb += X(i,j)*Beta[j];
    }
    TEMP += (Y[i] - xb) * X.row(i);
  }
  return(TEMP / n);
}
