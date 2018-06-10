#ifndef _SUBSET_H
#define _SUBSET_H

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix Subset(IntegerVector delta, NumericMatrix X, int n) {
	IntegerVector q=delta*seq_len(delta.size());
	IntegerVector set=q[q!=0];
	NumericMatrix nX(n,set.size());
	for(int i = 0; i < set.size(); i++) {
		nX(_,i)=X(_,set[i]-1);
	}
	return nX;
}

#else
#endif
