#include <RcppArmadillo.h>
#include "Basis.h"

using namespace Rcpp;

//' Generate a basis matrix
//'
//' @param X matrix.
//' @param knM lisy.
//' @examples
//' CombBasis(cbind(1:10,16:25),list(c(2,5,8),c(18,20,22)))

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
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
