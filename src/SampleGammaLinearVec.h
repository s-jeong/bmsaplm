#ifndef _SAMPLEGAMMALINEARVEC_H
#define _SAMPLEGAMMALINEARVEC_H

#include <RcppArmadillo.h>
#include "SampleGammalinear.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List SampleGammaLinearVec(NumericMatrix X, arma::vec y, int n, double a, List listdelta, IntegerVector gamma, int p, int q, double logpost, double dp){
	IntegerVector temp_gamma;
	double temp_logpost;
	temp_gamma = gamma;
	temp_logpost = logpost;
	List sam_gamma(2);
	for(int z = p; z < p+q; z++) {
		sam_gamma = SampleGammaLinear(X,y,n,a,listdelta,temp_gamma,p,q,z+1,temp_logpost,dp);
		temp_gamma = sam_gamma(0);
		temp_logpost = sam_gamma(1);
	}
	return sam_gamma;
}

#else
#endif
