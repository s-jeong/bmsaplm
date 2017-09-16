#ifndef _BASIS_H
#define _BASIS_H

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat Basis(arma::vec t, arma::vec knM) {
	arma::mat B = pow(abs(kron(t,arma::ones(1,knM.n_elem))-trans(kron(knM,arma::ones(1,t.n_elem)))),3);
	arma::mat jB = join_rows(t,B);
	return jB;
}

#else
#endif
