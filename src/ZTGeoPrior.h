#ifndef _ZTGeoPrior_H
#define _ZTGeoPrior_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double ZTGeoPrior(IntegerVector fordelta, int numkn, double phat) {
	IntegerVector seq_ind=seq_len(numkn);
	NumericVector pr = log(0.5)+dgeom(as<NumericVector>(wrap(sum(fordelta)-fordelta[0])), phat, true)-log(sum(dgeom(seq_ind, phat)))-lgamma(numkn+1)+lgamma(sum(fordelta)-fordelta[0]+1)+lgamma(numkn-(sum(fordelta)-fordelta[0])+1);
	return as<double>(wrap(pr));
}

// Be careful with log(0.5)

#else
#endif

