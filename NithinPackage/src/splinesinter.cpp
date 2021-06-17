// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
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

//' Check Spearman correlations between interactions in X and treatment
//' 
//' @param y A vector of outcomes.
//' @param X A matrix of spline bases.
//' @param alpha.schedule The prior on lambda
//' @export

//[[Rcpp::export]]
mat bs2(vec x, unsigned int deg) {
  BSpline b1{x, deg};
  BSpline b2{-x, deg};
  
  mat b = b1.basis();
  mat c = b2.basis();
  vec l = zeros(b.n_rows);
  int row_size = c.n_cols;
  for (int i = 0; i < b.n_rows; ++i) {
    l(i) = c(i, row_size - 1);
  }
  
  mat d = join_rows(l, b);
  return d;
}
//[[Rcpp::export]]
mat dbs2(vec x, unsigned int deg) {
  BSpline b1{x, deg};
  BSpline b2{-x, deg};
  mat db1 = b1.derivative();
  mat db2 = b2.derivative();
  vec l = zeros(db1.n_rows);
  int row_size = db2.n_cols;
  for (int i = 0; i < db1.n_rows; ++i) {
    l(i) = db2(i, row_size - 1);
  }
  
  mat d = join_rows(l, db1);
  return d;
}
//[[Rcpp::export]]
mat bsme(vec x) {
  
  mat m = join_rows(bs2(x, 3), bs2(x, 4), bs2(x, 5), bs2(x, 6));
  return m;
}
//[[Rcpp::export]]
mat dbsme(vec x) {
  
  mat m = join_rows(dbs2(x, 3), dbs2(x, 4), dbs2(x, 5), dbs2(x, 6));
  return m;
}

mat makebs(vec x) {
  vec v = ones(x.n_rows);
  return join_rows(v, bsme(x));
}

mat splineBases(vec x) {
  mat Xbs;
  for (int i = 0; i < X.n_rows; ++i) {
    Xbs = join_rows(Xbs, makebs(X));
  }
  return Xbs;
}

/*//[[Rcpp::export]]
void subSamp(vec v) {
  auto ve = RcppArmadillo::sample(v, v.size()/2, false);
  ve.print();
}*/

//[[Rcpp::export]]
struct Comp {
  
  public:
    bool operator()(vec a, vec b) {
      return a(3) < b(3);
    }
  
};

//[[Rcpp::export]]
List function(int obs, int covs, mat Xbs, vec y, vec treat, int k) {//k is number of top results, i.e. top 100 or top 300
  
  mat treatbs = bsme(treat);
  auto subsamp = sample(v, n/2, false);
  priority_queue<vec, vector<vec>, Comp> pq;
  
  for (int i = 0; i < treatbs.n_cols; ++i) {
    for (int j = 0; j < Xbs.n_cols; ++j) {
      for (int k = j; k < Xbs.n_cols; ++k) {
        double inter_temp = treatbs(subsamp, i)*Xbs(subsamp, j)*Xbs(subsamp, k);
        double cor_temp = cor(y[subsamp], inter_temp);
        vec indexCurr{i, j, k, cor_temp};
        pq.push(indexCurr);
        if (pq.size() > k) {
          pq.pop();
        }
      }
    }
  }
  List L;
  while (!pq.empty()) {
    L.push_back(pq.top());
    pq.pop();
  }
  return L;
}


//[[Rcpp::export]]
int main() {
  vec x{1,2,3,4};
  subSamp(x);
  return 0;
}

/*** R
main()
*/