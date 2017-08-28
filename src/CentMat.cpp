#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]

arma::mat CentMat(arma::mat X, arma::vec meanX) {
	arma::mat Y = X.each_row() - meanX.t();
	return Y;
}
