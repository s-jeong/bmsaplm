#ifndef _COMBBASIS_H
#define _COMBBASIS_H

#include <RcppArmadillo.h>
#include "Basis.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat CombBasis(arma::mat X, List knM) {
	arma::mat combB;
	for(int i = 0; i < as<int>(wrap(X.n_cols)); i++) {
		arma::vec t = X.col(i);
		arma::mat B = Basis(t,as<arma::vec>(knM[i]));
		combB = join_rows(combB,B);
	}
	return combB;
}

#else
#endif
