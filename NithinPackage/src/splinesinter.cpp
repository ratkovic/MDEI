// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppArmadilloExtensions/sample.h>
#include "splines2Armadillo.h"
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
  vec v = ones(x.size());
  mat m = join_rows(v, dbs2(x, 4), dbs2(x, 5), dbs2(x, 6));
  return m;
}
//[[Rcpp::export]]
mat makebs(mat X) {
  
  for (int i = 0; i < X.n_cols; ++i) {
    X.col(i) = (X.col(i) - mean(X.col(i)))/stddev(X.col(i));
  }
  
  vec x = X.as_col();
  vec v = ones(x.n_rows);
  map<int, int> m;
  for (int i = 0; i < x.size(); ++i) {
    m[x(i)]++;
  }
  if (m.size() <= 3) {
    return join_rows(v, x);
  }
  return join_rows(v, bsme(x));
}
//[[Rcpp::export]]
mat splineBases(mat X, int covs) {
  mat Xbs;
  for (int i = 0; i < covs; ++i) {
    Xbs = join_rows(Xbs, makebs(X));
  }
  return Xbs;
}

//[[Rcpp::export]]
vec subSamp(vec v) {
  return Rcpp::RcppArmadillo::sample(v, v.size()/2, false);
}

struct Comp { //this is a comparator, used for the heap (priority_queue) in the function below
  
  public:
    bool operator()(vec a, vec b) {
      return a(3) > b(3);
    }
  
};

//[[Rcpp::export]]
List correlations(int obs, int covs, mat X, vec y, vec treat, long long unsigned int a) {//a is number of top results, i.e. top 100 or top 300
  
  vec v = ones(obs);
  for (int i = 0; i < obs; ++i) {
    v(i) = i;
  }
  
  mat Xbs = splineBases(X, covs);
  
  mat treatbs = bsme(treat);
  vec sample = subSamp(v);
  
  mat treatSubsamp = zeros(obs/2, treatbs.n_cols);
  mat XSubsamp = zeros(obs/2, Xbs.n_cols);
  vec ySubsamp = zeros(obs/2);
  
  for (int i = 0; i < obs/2; ++i) { //this just gets the certain rows/entries corresponding to the random sample
    for (unsigned int j = 0; j < treatbs.n_cols; ++j) {
      treatSubsamp(i, j) = treatbs(sample(i), j);
      XSubsamp(i, j) = Xbs(sample(i), j);
      ySubsamp(i) = y(sample(i));
    }  
  }
  
  priority_queue<vec, vector<vec>, Comp> pq;
  double sdy = stddev(ySubsamp);

  for (double i = 0; i < treatbs.n_cols; ++i) {
    for (double j = 0; j < Xbs.n_cols; ++j) {
      for (double k = j; k < Xbs.n_cols; ++k) {
        vec inter_temp = treatSubsamp.col(i) % XSubsamp.col(j) % XSubsamp.col(k); 
        
        double cor_temp = 0;
        if (!inter_temp.is_zero()) {
          double cor_temp = abs(as_scalar(cor(ySubsamp, inter_temp)));
        }
        
        vec indexCurr{i, j, k, cor_temp};
        pq.push(indexCurr);
        if (pq.size() > a) { //this keeps the size of the heap to exactly a
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
  
  mat M;

  for (int i = 0; i < L.size(); ++i) {
    vec indexCurr = L[i];
    vec interTemp = treatSubsamp.col(indexCurr(0)) % XSubsamp.col(indexCurr(1)) % XSubsamp.col(indexCurr(2));
    M = join_rows(M, interTemp);
  }
  
  L.push_back(M);
  
  return L; //returns a highest correlations
} 

//[[Rcpp::export]]
vec checkcor(mat cors, double thresh) {
  vec v = ones(cors.n_cols); //include all vars initially. This is a bitmask where 1 means to include var and 0 means not to
  for (int i = 0; i < cors.n_rows; ++i) {
    for (int j = i; j < cors.n_cols; ++j) {
      if (abs(cors(i, j)) > thresh) {
        v(j) = 0;
      }
    }
  }
  return v; //vars marked zero are ones to not include
}
 
//[[Rcpp::export]]
mat gramschmidt(vec y, mat X) { //not finished yet. Right now it just finds the column with the max correlation
  double m = 0;
  double index = -1;
  for (int i = 0; i < X.n_cols; ++i) {
    vec col = X.col(i);
    if (as_scalar(cor(col, y)) > m) {
      index = i;
      m = as_scalar(cor(col, y));
    }
  }
}

//[[Rcpp::export]]
int main() {
  
  mat X = randn(1000, 5);
  vec treat = randn(1000);
  vec y = randn(1000);
  //bsme(y);
  clock_t a = clock();
  List L = correlations(1000, 5, X, y, treat, 100);
  clock_t b = clock();
  cout << double(b - a)/CLOCKS_PER_SEC << endl;
  return 0;
}

/*** R
main()
*/