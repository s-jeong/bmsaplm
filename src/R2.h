#ifndef _R2_H
#define _R2_H

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double R2(arma::mat X, arma::vec y) {
//	arma::vec bhat = solve(X.t()*X,X.t()*y);
	arma::vec bhat = inv_sympd(X.t()*X)*X.t()*y;
	double vR2 = sum(pow(X*bhat,2))/sum(pow(y-mean(y),2));
	return vR2;
}

#else
#endif
