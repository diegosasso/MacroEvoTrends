#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//' Matrix Exponential (scalar implementation)
//'
//' @param t time
//' @param A rate matrix
//'
//' @return probability matrix
//' @export
// [[Rcpp::export]]
arma::mat ExpMat_C(double t, arma::mat A) {
  arma::mat B = expmat(A*t);
  return(B);
}