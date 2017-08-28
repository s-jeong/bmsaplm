#ifndef _MAINSAMPLEDELTAGLOBALUPDATE_H
#define _MAINSAMPLEDELTAGLOBALUPDATE_H

#include <RcppArmadillo.h>
#include "MainSampleDeltaUpdate.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MainSampleDeltaGlobalUpdate(NumericMatrix X, arma::vec y, List indset, List delta, IntegerVector gamma, List listdelta, int numkn, int n, int p, int q, int z, double logbf_cur, double phat) {
	List delta_new = delta;
	double logbf_new = logbf_cur;
	for(int i = 0; i < indset.size(); i++) {
		IntegerVector sel_ind=indset(i);
		List list_samdelta = MainSampleDeltaUpdate(X,y,sel_ind,delta_new,gamma,listdelta,numkn,n,p,q,z,logbf_new,phat);
		delta_new = list_samdelta(0);
		logbf_new = as<double>(list_samdelta(1));
	}
	return List::create(delta_new,logbf_new);
}

#else
#endif

