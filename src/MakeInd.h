#ifndef _MAKEIND_H
#define _MAKEIND_H

#include <RcppArmadillo.h>
//#include <RcppArmadilloExtensions/sample.h>
//#include "RcppArmadilloSample.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MakeInd(int numkn) {
	IntegerVector ind;
	IntegerVector cs;
	do{
		ind = Rcpp::RcppArmadillo::sample(IntegerVector::create(2,3,4),numkn,TRUE);
//	  ind = Rcpp::RcppArmadillo::sample(IntegerVector::create(2,3,4),numkn,TRUE,NumericVector::create(1/3,1/3,1/3));
		IntegerVector sugar_cs = cumsum(ind);
		cs = sugar_cs;
	}while(is_false(any(cs==numkn)));
	IntegerVector sq = seq_len(numkn);
	int k = as<int>(wrap(sq[cs==numkn]));
	List ind_set(k+1);
	ind_set(0)=1;
	for(int i = 0; i < k; i++) {
		ind_set(i+1)=Range(cs[i]-ind[i]+2,cs[i]+1);
	}
	return ind_set;
}

#else
#endif
