// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
#include "splines2Armadillo.h"
#include <queue>
#include <vector>

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
  for (unsigned int i = 0; i < b.n_rows; ++i) {
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
  for (unsigned int i = 0; i < db1.n_rows; ++i) {
    l(i) = db2(i, row_size - 1);
  }
  
  mat d = join_rows(l, db1);
  return d;
}
//[[Rcpp::export]]
mat bsme(vec x) {
  
  mat m = join_rows(x, bs2(x, 4), bs2(x, 5), bs2(x, 6));
  return m;
}
//[[Rcpp::export]]
mat dbsme(vec x) {
  
  mat m = join_rows(x, dbs2(x, 4), dbs2(x, 5), dbs2(x, 6));
  return m;
}

mat makebs(vec x) {
  vec v = ones(x.n_rows);
  return join_rows(v, bsme(x));
}

mat splineBases(vec x, int covs) {
  mat Xbs;
  for (int i = 0; i < covs; ++i) {
    Xbs = join_rows(Xbs, makebs(x));
  }
  return Xbs;
}

//[[Rcpp::export]]
vec subSamp(vec v) {
  return Rcpp::RcppArmadillo::sample(v, v.size()/2, false);
}

struct Comp {
  
  public:
    bool operator()(vec a, vec b) {
      return a(3) < b(3);
    }
  
};

//[[Rcpp::export]]
List correlations(int obs, int covs, vec x, vec y, vec treat, long long unsigned int a) {//a is number of top results, i.e. top 100 or top 300
  
  vec v = ones(obs);
  for (int i = 0; i < obs; ++i) {
    v(i) = i;
  }
  
  mat Xbs = splineBases(x, covs);
  
  mat treatbs = bsme(treat);
  vec sample = subSamp(v);
  
  mat treatSubsamp = zeros(obs/2, treatbs.n_cols);
  mat XSubsamp = zeros(obs/2, Xbs.n_cols);
  vec ySubsamp = zeros(obs/2);
  
  for (int i = 0; i < obs/2; ++i) {
    for (unsigned int j = 0; j < treatbs.n_cols; ++j) {
      treatSubsamp(i, j) = treatbs(sample(i), j);
      XSubsamp(i, j) = Xbs(sample(i), j);
      ySubsamp(i) = y(sample(i));
    }  
  }
  
  //treatSubsamp = treatSubsamp.t();
  //XSubsamp = XSubsamp.t();
  
  //cout << treatSubsamp(0, 0) << endl;
  //cout << XSubsamp(0) << endl;
  
  priority_queue<vec, vector<vec>, Comp> pq;
  vec zero_vec = zeros(Xbs.n_cols);
  
  for (double i = 0; i < treatbs.n_cols; ++i) {
    for (double j = 0; j < Xbs.n_cols; ++j) {
      for (double k = j; k < Xbs.n_cols; ++k) {
        vec inter_temp = treatSubsamp.col(i) % XSubsamp.col(j) % XSubsamp.col(k); 
        double cor_temp = 0;
        if (approx_equal(zero_vec, inter_temp, "absdiff", 0.001)) {
          cor_temp = as_scalar(cor(ySubsamp, inter_temp));
        }
        vec indexCurr{i, j, k, cor_temp};
        pq.push(indexCurr);
        if (pq.size() > a) {
          pq.pop();
        }
      }
    }
  }
  List L;
  while (!pq.empty()) {
    L.push_back(pq.top());
    pq.pop();
    cout << pq.top()(3) << endl;
  }
  return L; //returns k highest correlations
} 


//[[Rcpp::export]]
int main() {
  vec x{3,1,4,1};
  correlations(4, 2, x, x, x, 2);
  return 0;
}

/*** R
main()
*/