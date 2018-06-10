#ifndef _ETA_H
#define _ETA_H

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

IntegerVector Eta(List delta, IntegerVector gamma, int p, int q) {
	int total_length = 0;
	for (int i = 0; i < p; ++i){
		total_length += Rf_length(delta[i]);
	}
	IntegerVector ind(p);
	ind[gamma[Range(0,p-1)]==2] = 1;
	int index = 0;
	IntegerVector vecdelta(total_length+q);
	for(int i = 0; i < p; i++) {
		IntegerVector tempdelta=as<IntegerVector>(wrap(delta[i]))*ind(i);
		if(gamma[i]==1) {
			tempdelta[0]=1;
		}
		std::copy(tempdelta.begin(), tempdelta.end(), vecdelta.begin() + index);
		index += tempdelta.size();
	}
	if(q > 0.5){
	  std::copy(gamma[Range(p,p+q-1)].begin(), gamma[Range(p,p+q-1)].end(), vecdelta.begin() + index);
	}
	return vecdelta;
}

#else
#endif
