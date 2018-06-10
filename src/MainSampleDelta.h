#ifndef _MAINSAMPLEDELTA_H
#define _MAINSAMPLEDELTA_H

#include <RcppArmadillo.h>
// #include <RcppArmadilloExtensions/sample.h>
// #include "RcppArmadilloSample.h"
#include "ZTGeoPrior.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

IntegerVector MainSampleDelta(IntegerVector sel_ind, IntegerVector fordelta, List listdelta, int numkn, double phat) {
	NumericVector prob(pow(2,sel_ind.size()));
	IntegerVector curdelta(numkn+1);
	for(int i = 0; i < numkn+1; i++) {
		curdelta[i]=fordelta[i];
	}
	List templistdelta=listdelta[sel_ind.size()-1];
	for(int i = 0; i < pow(2,sel_ind.size()); i++) {
		NumericVector cand = templistdelta[i];
		for(int j = 0; j < sel_ind.size(); j++) {
			curdelta[sel_ind[j]-1] = cand[j];
		}
//		prob[i]=ZTBinPrior(curdelta,numkn,phat);
		prob[i]=ZTGeoPrior(curdelta,numkn,phat);
	}
	CharacterVector sel_ch=as<CharacterVector>(wrap(seq_len(pow(2,sel_ind.size()))));
	IntegerVector seq_ind=seq_len(pow(2,sel_ind.size()));
	do{
		IntegerVector sampled_ind = Rcpp::RcppArmadillo::sample(seq_ind,1,TRUE,exp(prob-max(prob)));
//		IntegerVector sampled_ind = RcppArmadilloSample(seq_ind,1,TRUE,exp(prob-max(prob)));
		IntegerVector sampled_val = templistdelta[sampled_ind[0]-1];
		for(int j = 0; j < sel_ind.size(); j++) {
			curdelta[sel_ind[j]-1] = sampled_val[j];
		}
	}while(sum(curdelta)-curdelta[0]==0);
	return curdelta;
}

#else
#endif
