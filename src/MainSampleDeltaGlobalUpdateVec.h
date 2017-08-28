#ifndef _MAINSAMPLEDELTAGLOBALUPDATEVEC_H
#define _MAINSAMPLEDELTAGLOBALUPDATEVEC_H

#include <RcppArmadillo.h>
#include "MainSampleDeltaGlobalUpdate.h"
#include "MakeInd.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MainSampleDeltaGlobalUpdateVec(NumericMatrix X, arma::vec y, List delta, IntegerVector gamma, List listdelta, IntegerVector numkn, int n, int p, int q, double logbf_cur, NumericVector phat) {
	List temp_delta(p);
	double temp_logbf;
	temp_delta = delta;
	temp_logbf = logbf_cur;
	List sam_delta(2);
	for(int z = 0; z < p; z++) {
		if(gamma[z]==2){
			List indset = MakeInd(numkn[z]);
			sam_delta = MainSampleDeltaGlobalUpdate(X,y,indset,temp_delta,gamma,listdelta,numkn[z],n,p,q,z+1,temp_logbf,phat[z]);
			temp_delta = sam_delta(0);
			temp_logbf = sam_delta(1);
		}
	}
	return sam_delta;
}

#else
#endif
