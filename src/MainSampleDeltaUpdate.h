#ifndef _MAINSAMPLEDELTAUPDATE_H
#define _MAINSAMPLEDELTAUPDATE_H

#include <RcppArmadillo.h>
#include "MainSampleDelta.h"
#include "Eta.h"
#include "LogBF.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MainSampleDeltaUpdate(NumericMatrix X, arma::vec y, IntegerVector sel_ind, List delta, IntegerVector gamma, List listdelta, int numkn, int n, int p, int q, int z, double logbf_cur, double phat) {
	IntegerVector fordelta=delta[z-1];
	IntegerVector fordelta_cand=MainSampleDelta(sel_ind,fordelta,listdelta,numkn,phat);
	List delta_cand(p);
	List delta_new(p);
	double logbf_new;
	if(is_false(all(fordelta_cand==fordelta))){
		for(int i = 0; i < p; i++) {
			if(i==z-1){
				delta_cand[i]=fordelta_cand;
			}else{
				delta_cand[i]=delta[i];
			}
		}
		IntegerVector vdelta=Eta(delta_cand,gamma,p,q);
		double logbf_cand=LogBF(vdelta,X,y,n,-0.75);
		if(as<double>(runif(1))<exp(logbf_cand-logbf_cur)){
			delta_new=delta_cand;
			logbf_new=logbf_cand;
		}else{
			delta_new=delta;
			logbf_new=logbf_cur;
		}
	}else{
		delta_new=delta;
		logbf_new=logbf_cur;
	}
	return List::create(delta_new,logbf_new);
}

#else
#endif
