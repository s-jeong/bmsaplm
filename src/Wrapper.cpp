#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "CombBasis.h"
#include "CentMat.h"
#include "MCMCEstIteration.h"
#include "MCMCOneIteration.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]

arma::mat CombBasisCpp(arma::mat X, List knM) {
  arma::mat combB = CombBasis(X, knM);
  return combB;
}

arma::mat CentMatCpp(arma::mat X, arma::vec meanX) {
  arma::mat Y = CentMat(X, meanX);
  return Y;
}

List MCMCEstIterationCpp(NumericMatrix X, arma::vec y, List delta, IntegerVector gamma, List listdelta, IntegerVector numkn, int n, int p, int q, double logbf_cur,double a) {
  List output1 = MCMCEstIteration(X, y, delta, gamma, listdelta, numkn, n, p, q, logbf_cur,a);
  return output1;
}

List MCMCOneIterationCpp(NumericMatrix X, arma::vec y, List delta, IntegerVector gamma, List listdelta, IntegerVector numkn, int n, int p, int q, double logbf_cur, NumericVector phat, NumericVector mu, NumericVector var, NumericVector mu0, double a) {
  List output2 = MCMCOneIteration(X, y, delta, gamma, listdelta, numkn, n, p, q, logbf_cur, phat, mu, var, mu0, a);
  return output2;
}
