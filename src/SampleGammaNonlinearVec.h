#ifndef _SAMPLEGAMMANONLINEARVEC_H
#define _SAMPLEGAMMANONLINEARVEC_H

#include <RcppArmadillo.h>
#include "SampleGammaNonlinear.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List SampleGammaNonlinearVec(NumericMatrix X, arma::vec y, int n, double a, List listdelta, IntegerVector gamma, NumericVector mu, NumericVector var, IntegerVector numkn, NumericVector mu0, NumericVector phat, int p, int q, double logpost) {
	IntegerVector temp_gamma;
	double temp_logpost;
	temp_gamma = gamma;
	temp_logpost = logpost;
	List sam_gamma(2);
	for(int z = 0; z < p; z++) {
		sam_gamma = SampleGammaNonlinear(X,y,n,a,listdelta,temp_gamma,mu,var,numkn,mu0,phat,p,q,z+1,temp_logpost);
		temp_gamma = sam_gamma(0);
		temp_logpost = sam_gamma(1);
	}
	return sam_gamma;
}

#else
#endif
