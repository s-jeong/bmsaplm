#ifndef _LOGBF_H
#define _LOGBF_H

#include <RcppArmadillo.h>
#include "R2.h"
#include "Subset.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double LogBF(IntegerVector delta, NumericMatrix X, arma::vec y, int n, double a) {
	int q = sum(delta);
	double LogBF;
	if(q==0){
		LogBF = 0;
	}else{
		NumericMatrix dX = Subset(delta,X,n);
		double R2v = R2(as<arma::mat>(dX),y);
		LogBF = Rf_lbeta(0.5*q+a+1,0.5*(n-q-3)-a)-Rf_lbeta(a+1,0.5*(n-q-3)-a)-(0.5*(n-q-3)-a)*log(1-R2v);
	}
	return LogBF;
}

#else
#endif
