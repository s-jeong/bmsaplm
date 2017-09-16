#ifndef _SAMPLEFIXEFF_H
#define _SAMPLEFIXEFF_H

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat SampleFixeff(int n, arma::mat WXi, arma::vec y, double qv){
	arma::mat XtX = WXi.t()*WXi;
	arma::mat inv_XtX = inv_sympd(XtX);
	arma::mat Sigma = qv/(qv+1)*inv_XtX;
	arma::mat mu = qv/(qv+1)*inv_XtX*WXi.t()*(y-mean(y));
	int p = Sigma.n_cols;
	arma::mat X = reshape(arma::vec(rnorm(p*n)),p,n);
	arma::vec eigval;
	arma::mat eigvec;
	eig_sym(eigval,eigvec,Sigma);
	X = eigvec * diagmat(sqrt(eigval))*X;
	X.each_col() += mu;
	return X.t();
}

#else
#endif
