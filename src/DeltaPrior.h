#ifndef _DELTAPRIOR_H
#define _DELTAPRIOR_H

#include <RcppArmadillo.h>
#include "PseudoPrior.h"
#include "ZTGeoPrior.h"


using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double DeltaPrior(List listdelta, IntegerVector gamma, NumericVector mu, NumericVector var, IntegerVector numkn, NumericVector mu0, NumericVector phat, int p) {
	NumericVector pr(p);
	for(int i = 0; i < p; i++) {
		if(gamma[i]==2){
//			pr[i]=ZTBinPrior(listdelta[i],numkn[i],phat[i]);
			pr[i]=ZTGeoPrior(listdelta[i],numkn[i],phat[i]);
		}else{
			pr[i]=PseudoPrior(listdelta[i],mu[i],var[i],numkn[i],mu0[i]);
		}
	}
	double sumpr=sum(pr);
	return sumpr;
}

#else
#endif
