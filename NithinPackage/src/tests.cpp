// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
#include "splines2Armadillo.h"
#include <splinesinter.cpp>
#include <queue>
#include <vector>
#include <time.h>
#include <map>

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace splines2;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(splines2)]]

//' Check Spearman correlations between interactions in X and treatment
//' 
//' @param y A vector of outcomes.
//' @param X A matrix of spline bases.
//' @param alpha.schedule The prior on lambda
//' @export



//[[Rcpp::export]]
int main() {
  
  //tests
  vec a = randn(10);
  arma::mat b = bs2(a, 4);
  b.print();
  cout << "\n" << endl;
  arma::mat c = dbs2(a, 5);
  c.print();
  
  cout << "\n" << endl;
  
  arma::mat d = bsme(a);
  d.print();
  
  cout << "\n" << endl;
  
  arma::mat e = dbsme(a);
  e.print();
  
  cout << "\n" << endl;
  
  arma::mat f = randn(5,5);
  arma::mat g = makebs(f);
  
  g.print();
  
  cout << "\n" << endl;
  arma::mat h = splineBases(f, 5);
  h.print(); 

  /*arma::mat X = randn(1000, 5);
  //splineBases(X, 2);
  arma::vec treat = randn(1000);
  arma::vec y = randn(1000);
  clock_t a = clock();
  List L = correlations(1000, 5, X, {"education", "race", "ethnicity", "income", "other"}, y, treat, "treatname", 10);
  clock_t b = clock();
  cout << double(b - a)/CLOCKS_PER_SEC << endl; */
  return 0;
}

/*** R
main()
*/