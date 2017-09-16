#ifndef _ZTBINPRIOR_H
#define _ZTBINPRIOR_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double ZTBinPrior(IntegerVector fordelta, int numkn, double phat) {
	double pr = log(0.5)+(sum(fordelta)-fordelta[0])*log(phat)+(numkn-(sum(fordelta)-fordelta[0]))*log(1-phat)-log(1-pow(1-phat,numkn));
	return pr;
}

// Be careful with log(0.5)

#else
#endif

