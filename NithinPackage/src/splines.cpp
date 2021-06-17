// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "splines2Armadillo.h"

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace splines2;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(splines2)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//


//' Check Spearman correlations between interactions in X and treatment
//' 
//' @param y A vector of outcomes.
//' @param X A matrix of spline bases.
//' @param alpha.schedule The prior on lambda
//' @export
//[[Rcpp::export]]

void splineMatrix(vec y) {
  
  vec x = ones(5); //initializing to size 100
  for (int i = 0; i < 5; ++i) {
    x(i) = i/5.0;
  }
  BSpline splineMat{x, 5};
  mat b = splineMat.basis();
  mat c = cor(b, y);
  c.print();
  
}
//[[Rcpp::export]]
int fib(int n) {
  if (n < 2) {
    return n;
  }
  return fib(n - 1) + fib(n - 2);
}

//[[Rcpp::export]]
int main() {
  vec y{3,1,4,2,5};
  splineMatrix(y);
  //cout << fib(12) << endl;
  return 0;
}

/*** R
main()
*/