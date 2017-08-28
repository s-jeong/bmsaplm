#ifndef _ZTBBPRIOR_H
#define _ZTBBPRIOR_H

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double ZTBBPrior(IntegerVector fordelta, int numkn) {
	double pr = Rf_lbeta(sum(fordelta)-fordelta[0]+1,numkn-(sum(fordelta)-fordelta[0])+1);
	return pr;
}

#else
#endif
