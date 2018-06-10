#ifndef _SAMPLEGAMMANONLINEAR_H
#define _SAMPLEGAMMANONLINEAR_H

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

List SampleGammaNonlinear(NumericMatrix X, arma::vec y, int n, double a, List listdelta, IntegerVector gamma, NumericVector mu, NumericVector var, IntegerVector numkn, NumericVector mu0, NumericVector phat, int p, int q, int z, double logpost) {
	NumericVector vec_logpost(3);
	IntegerVector temp_gamma(p+q);
	temp_gamma[Range(0,p+q-1)]=gamma;
	for(int i = 0; i < 3; i++) {
		if(gamma[z-1]==i){
			vec_logpost[i] = logpost;
		}else{
			temp_gamma[z-1] = i;
			double gp = GammaPrior(temp_gamma,p,q);
			double dp = DeltaPrior(listdelta,temp_gamma,mu,var,numkn,mu0,phat,p);
			IntegerVector vdelta = Eta(listdelta,temp_gamma,p,q);
			double lbf = LogBF(vdelta,X,y,n,a);
			vec_logpost[i] = lbf+gp+dp;
		}
	}
	IntegerVector seq_ind = seq_len(3)-1;
	IntegerVector sampled_gamma = Rcpp::RcppArmadillo::sample(seq_ind,1,TRUE,exp(vec_logpost-max(vec_logpost)));
	temp_gamma[z-1] = sampled_gamma[0];
	return List::create(temp_gamma,vec_logpost[sampled_gamma[0]]);
}

#else
#endif
