#ifndef _SAMPLEGAMMALINEAR_H
#define _SAMPLEGAMMALINEAR_H

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "GammaPrior.h"
#include "DeltaPrior.h"
#include "Eta.h"
#include "LogBF.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List SampleGammaLinear(NumericMatrix X, arma::vec y, int n, double a, List listdelta, IntegerVector gamma, int p, int q, int z, double logpost, double dp) {
	NumericVector vec_logpost(2);
	IntegerVector temp_gamma(p+q);
	temp_gamma[Range(0,p+q-1)]=gamma;
	for(int i = 0; i < 2; i++) {
		if(gamma[z-1]==i){
			vec_logpost[i] = logpost;
		}else{
			temp_gamma[z-1] = i;
			double gp = GammaPrior(temp_gamma,p,q);
			IntegerVector vdelta = Eta(listdelta,temp_gamma,p,q);
			double lbf = LogBF(vdelta,X,y,n,a);
			vec_logpost[i] = lbf+gp+dp;
		}
	}
	IntegerVector seq_ind = seq_len(2)-1;
	IntegerVector sampled_gamma = Rcpp::RcppArmadillo::sample(seq_ind,1,TRUE,exp(vec_logpost-max(vec_logpost)));
	temp_gamma[z-1] = sampled_gamma[0];
	return List::create(temp_gamma,vec_logpost[sampled_gamma[0]]);
}

#else
#endif
