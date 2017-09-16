#ifndef _PSEUDOSAMPLEDELTA_H
#define _PSEUDOSAMPLEDELTA_H

#include <RcppArmadillo.h>
// #include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

IntegerVector PseudoSampleDelta(double mu, double var, int numkn, double mu0) {
	IntegerVector seq_ind=seq_len(numkn);

	NumericVector prob = dnorm(seq_ind, mu, pow(var,0.5), true);
	IntegerVector sampled_ind = Rcpp::RcppArmadillo::sample(seq_ind,1,TRUE,exp(prob-max(prob)));
	IntegerVector sampled_loc = Rcpp::RcppArmadillo::sample(seq_ind,as<int>(wrap(sampled_ind)),FALSE);
	IntegerVector delta(numkn+1);
	for(int j = 0; j < sampled_loc.size(); j++) {
		delta[sampled_loc[j]] = 1;
	}
	if(as<double>(runif(1))<mu0){
		delta[0]=1;
	}

	return delta;
}

#else
#endif
