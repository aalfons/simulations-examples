// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat tcov_cpp(const arma::mat& x, const double& beta) {

  // initializations
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;

  // In the paper, we have w(x) = exp(-x/2). But since we always call
  // w(beta * r^2), we instead set b = -beta/2 and use w(x) = exp(x).
  const double b = -beta / 2.0;

  // compute inverse of the sample covariance matrix
  const arma::mat cov_inv = arma::solve(arma::cov(x),
                                        arma::eye<arma::mat>(p, p));

  // loop over pairs of observations
  arma::uword i, j, k, l;             // running indices
  arma::vec diff(p);                  // difference of pair of observations
  arma::mat V(p, p, fill::zeros);     // scatter matrix
  double r_sq, w, denominator = 0.0;  // squared distance, weight, denominator
  for(i = 1; i < n; i++) {
    for(j = 0; j < i; j++) {
      // compute difference of current pair of observations
      for(k = 0; k < p; k++) diff(k) = x(i,k) - x(j,k);
      // compute squared pairwise Mahalanobis distance
      r_sq = 0.0;
      for(k = 0; k < p; k++) {
        for(l = 0; l < p; l++) {
          r_sq += diff(k) * cov_inv(k,l) * diff(l);
        }
      }
      // compute weight for current pair of observations
      w = exp(b * r_sq);
      // add weighted contribution of current pair of observations
      for(k = 0; k < p; k++) {
        // update diagonal elements
        V(k,k) += w * diff(k) * diff(k);
        // update offdiagonal elements
        for(l = 0; l < k; l++) {
          V(k,l) += w * diff(k) * diff(l);
          V(l,k) = V(k,l);
        }
      }
      // add weight of current pair of observations to denominator
      denominator += w;
    }
  }

  // return scatter matrix
  return V / denominator;
}
