#ifndef _PILOTSAMPLEDELTAGLOBALUPDATEVEC_H
#define _PILOTSAMPLEDELTAGLOBALUPDATEVEC_H

#include <RcppArmadillo.h>
#include "PilotSampleDeltaGlobalUpdate.h"
#include "MakeInd.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List PilotSampleDeltaGlobalUpdateVec(NumericMatrix X, arma::vec y, List delta, IntegerVector gamma, List listdelta, IntegerVector numkn, int n, int p, int q, double logbf_cur) {
	List temp_delta(p);
	double temp_logbf;
	temp_delta = delta;
	temp_logbf = logbf_cur;
	List sam_delta(2);
	for(int z = 0; z < p; z++) {
		if(gamma[z]==2){
			List indset = MakeInd(numkn[z]);
			sam_delta = PilotSampleDeltaGlobalUpdate(X,y,indset,temp_delta,gamma,listdelta,numkn[z],n,p,q,z+1,temp_logbf);
			temp_delta = sam_delta(0);
			temp_logbf = sam_delta(1);
		}
	}
	return sam_delta;
}

#else
#endif
