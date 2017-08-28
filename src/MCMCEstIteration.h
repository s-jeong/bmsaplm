#ifndef _MCMCESTITERATION_H
#define _MCMCESTITERATION_H

#include <RcppArmadillo.h>
// #include <RcppArmadilloExtensions/sample.h>
#include "PilotSampleDeltaGlobalUpdateVec.h"
#include "SampleFixeff.h"
#include "Eta.h"
#include "R2.h"
#include "SampleFixeff.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]

List MCMCEstIteration(NumericMatrix X, arma::vec y, List delta, IntegerVector gamma, List listdelta, IntegerVector numkn, int n, int p, int q, double logbf_cur,double a) {
	List sam_delta = PilotSampleDeltaGlobalUpdateVec(X,y,delta,gamma,listdelta,numkn,n,p,q,logbf_cur);
	List temp_delta = sam_delta(0);
	double temp_logbf = sam_delta(1);
	IntegerVector vdelta = Eta(temp_delta,gamma,p,q);
	arma::mat dX = as<arma::mat>(Subset(vdelta,X,n));
	double vR2 = R2(dX,y);
	double g_star = R::rbeta(0.5*sum(vdelta)+a+1,0.5*(n-sum(vdelta)-3)-a);
	double g = (1/g_star-1)/(1-vR2);
	double sigsq = 1/R::rgamma(0.5*(n-1),(2*(1+g))/(sum(pow(y-mean(y),2))*(1+g*(1-vR2))));	// note scale = 1/rate
	arma::mat fixeff = SampleFixeff(1,dX,y,g);
	double alpha0=R::rnorm(mean(y),pow(sigsq/n,0.5));
	return List::create(temp_delta,vdelta,temp_logbf,g,sigsq,fixeff,alpha0);
}

#else
#endif
