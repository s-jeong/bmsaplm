#ifndef _PSEUDOSAMPLEDELTAVEC_H
#define _PSEUDOSAMPLEDELTAVEC_H

#include <RcppArmadillo.h>
#include "PseudoSampleDelta.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List PseudoSampleDeltaVec(List delta, IntegerVector gamma, NumericVector mu, NumericVector var, IntegerVector numkn, NumericVector mu0, int p) {
	IntegerVector gamma2 = pmax(gamma-1,0);
	int k=sum(gamma2);
	List temp_delta(p-k);
	int m=0;
	for(int z = 0; z < p; z++) {
		if(gamma[z]!=2){
			temp_delta[m]=PseudoSampleDelta(mu[z],var[z],numkn[z],mu0[z]);
			m++;
		}
	}
	return temp_delta;
}

#else
#endif
