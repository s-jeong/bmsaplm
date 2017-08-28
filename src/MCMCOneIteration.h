#ifndef _MCMCONEITERATION_H
#define _MCMCONEITERATION_H

#include <RcppArmadillo.h>
#include "MainSampleDeltaGlobalUpdateVec.h"
#include "PseudoSampleDeltaVec.h"
#include "GammaPrior.h"
#include "DeltaPrior.h"
#include "SampleGammaNonlinearVec.h"
#include "SampleGammaLinearVec.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MCMCOneIteration(NumericMatrix X, arma::vec y, List delta, IntegerVector gamma, List listdelta, IntegerVector numkn, int n, int p, int q, double logbf_cur, NumericVector phat, NumericVector mu, NumericVector var, NumericVector mu0, double a) {
	List temp_delta(p);
	IntegerVector temp_gamma;
	double temp_logbf, temp_logpost;
	temp_delta = delta;
	temp_gamma = gamma;
	temp_logbf = logbf_cur;
	if(is_true(any(gamma==2))){
		List delta_iter=MainSampleDeltaGlobalUpdateVec(X,y,temp_delta,temp_gamma,listdelta,numkn,n,p,q,temp_logbf,phat);
		temp_delta = delta_iter(0);
	}
	if(is_true(any(gamma!=2))){
		List pse_delta = PseudoSampleDeltaVec(temp_delta,temp_gamma,mu,var,numkn,mu0,p);
		int m = 0;
		for(int z = 0; z < p; z++) {
			if(gamma[z]!=2){
				temp_delta[z] = pse_delta[m];
				m++;
			}
		}
	}
	IntegerVector vdelta=Eta(temp_delta,temp_gamma,p,q);
	double lprior_gamma=GammaPrior(temp_gamma,p,q);
	double lprior_delta=DeltaPrior(temp_delta,temp_gamma,mu,var,numkn,mu0,phat,p);
	temp_logbf = LogBF(vdelta,X,y,n,a);
	temp_logpost = temp_logbf+lprior_gamma+lprior_delta;
	List sam_gammaN = SampleGammaNonlinearVec(X,y,n,a,temp_delta,temp_gamma,mu,var,numkn,mu0,phat,p,q,temp_logpost);
	temp_gamma = sam_gammaN(0);
	temp_logpost = sam_gammaN(1);
	lprior_delta=DeltaPrior(temp_delta,temp_gamma,mu,var,numkn,mu0,phat,p);
	List sam_gammaL = SampleGammaLinearVec(X,y,n,a,temp_delta,temp_gamma,p,q,temp_logpost,lprior_delta);
	temp_gamma = sam_gammaL(0);
	temp_logpost = sam_gammaL(1);
	temp_logbf = temp_logpost-lprior_delta-GammaPrior(temp_gamma,p,q);
	return List::create(temp_delta,temp_gamma,temp_logbf);
}

#else
#endif

