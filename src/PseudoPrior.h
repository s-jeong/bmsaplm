#ifndef _PSEUDOPRIOR_H
#define _PSEUDOPRIOR_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double PseudoPrior(IntegerVector delta, double mu, double var, int numkn, double mu0) {
	IntegerVector seq_ind=seq_len(numkn);
	NumericVector pr = delta[0]*log(mu0)+(1-delta[0])*log(1-mu0)+dnorm(as<NumericVector>(wrap(sum(delta)-delta[0])), mu, pow(var,0.5), true)-log(sum(dnorm(seq_ind, mu, pow(var,0.5))))-lgamma(numkn+1)+lgamma(sum(delta)-delta[0]+1)+lgamma(numkn-(sum(delta)-delta[0])+1);
	return as<double>(wrap(pr));
}

#else
#endif
