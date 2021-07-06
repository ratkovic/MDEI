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
arma::mat bs2(arma::vec x, unsigned int deg) {
  BSpline b1{x, deg + 1};
  BSpline b2{-x, deg + 1};
  
  arma::mat b = b1.basis();
  arma::mat c = b2.basis();
  arma::vec l = arma::zeros(b.n_rows);
  int row_size = c.n_cols;
  for (unsigned int i = 0; i < b.n_rows; ++i) {
    l(i) = c(i, row_size - 1);
  }
  
  arma::mat d = arma::join_rows(l, b);
  return b;
}
//[[Rcpp::export]]
arma::mat dbs2(arma::vec x, unsigned int deg) {
  BSpline b1{x, deg + 1};
  BSpline b2{-x, deg + 1};
  arma::mat db1 = b1.derivative();
  arma::mat db2 = b2.derivative();
  arma::vec l = arma::zeros(db1.n_rows);
  int row_size = db2.n_cols;
  for (unsigned int i = 0; i < db1.n_rows; ++i) {
    l(i) = -db2(i, row_size - 1);
  }
  
  arma::mat d = arma::join_rows(l, db1);
  return db1;
}
//[[Rcpp::export]]
arma::mat bsme(arma::vec x) {
  
  arma::mat m = arma::join_rows(x, bs2(x, 3), bs2(x, 4), bs2(x, 5));
  return arma::join_rows(m, bs2(x, 6));
}
//[[Rcpp::export]]
arma::mat dbsme(arma::vec x) {
  arma::vec v = arma::ones(x.n_rows);
  arma::mat m = arma::join_rows(v, dbs2(x, 4), dbs2(x, 4), dbs2(x, 5));
  return arma::join_rows(m, dbs2(x, 6));
}
//[[Rcpp::export]]
arma::mat makebs(arma::vec X) {

  arma::vec x = (X - mean(X))/stddev(X);
  arma::vec v = arma::ones(x.n_rows);
  map<double, int> m;
  for (unsigned int i = 0; i < x.n_rows; ++i) {
    m[x(i)]++;
  }
  if (m.size() <= 3) {
    return arma::join_rows(v, x);
  }
  return arma::join_rows(v, bsme(x));
}
//[[Rcpp::export]]
List splineBases(arma::mat X, int covs) {
  arma::mat Xbs;
  std::vector<int> v;
  v.push_back(0);
  for (int i = 0; i < covs; ++i) {
    arma::mat m = makebs(X.col(i));
    Xbs = arma::join_rows(Xbs, m);
    v.push_back(m.n_cols + v[i]);

  }
  return List::create(Named("Xbs") = Xbs, _["vec"] = v);
}

//[[Rcpp::export]]
arma::vec subSamp(arma::vec v) {
  return Rcpp::RcppArmadillo::sample(v, v.size()/2, false);
}

struct Comp { //this is a comparator, used for the heap (priority_queue) in the function below
  public:
    bool operator()(arma::vec a, arma::vec b) {
      return a(3) > b(3);
    }
};

//[[Rcpp::export]]
List splineBasesAndCorrs(int obs, int covs, arma::mat X, std::vector<std::string> Xname, arma::vec y, arma::vec treat, std::string treatName, long long unsigned int a) {
  //a is number of top results, i.e. top 100 or top 300
  arma::vec v = arma::ones(obs);
  for (int i = 0; i < obs; ++i) {
    v(i) = i;
  }

  List sB = splineBases(X, covs);
  arma::mat Xbs = sB["Xbs"];
  std::vector<int> colSizes = sB["vec"];

  arma::mat treatbs = bsme(treat);
  arma::vec sample = subSamp(v);
  
  arma::mat treatSubsamp = arma::zeros(obs/2, treatbs.n_cols);
  arma::mat XSubsamp = arma::zeros(obs/2, Xbs.n_cols);
  arma::vec ySubsamp = arma::zeros(obs/2);
  
  for (int i = 0; i < obs/2; ++i) { //this just gets the certain rows/entries corresponding to the random sample
    for (unsigned int j = 0; j < treatbs.n_cols; ++j) {
      treatSubsamp(i, j) = treatbs(sample(i), j);
      XSubsamp(i, j) = Xbs(sample(i), j);
      ySubsamp(i) = y(sample(i));
    }  
  }

  priority_queue<arma::vec, std::vector<arma::vec>, Comp> pq;
  arma::vec indexCurr=arma::zeros(4);
  arma::vec inter_temp = arma::zeros(treatbs.n_rows);

  for (double i = 0; i < treatbs.n_cols; ++i) {
    for (double j = 0; j < Xbs.n_cols; ++j) {
      for (double k = j + 1; k < Xbs.n_cols; ++k) {
        inter_temp = treatSubsamp.col(i) % XSubsamp.col(j) % XSubsamp.col(k); 
        
        double cor_temp = 0;
        if (!inter_temp.is_zero()) {
          cor_temp = abs(as_scalar(arma::cor(ySubsamp, inter_temp)));
        }
        indexCurr(0) = i;
        indexCurr(1) = j;
        indexCurr(2) = k;
        indexCurr(3) = cor_temp;
        pq.push(indexCurr);
        if (pq.size() > a) { //this keeps the size of the heap to exactly a
          pq.pop();
        }
      }
    }
  }
  
  std::vector<arma::vec> indexCurrs;

  while (!pq.empty()) {
    indexCurrs.push_back(pq.top());
    pq.pop();
  }
  
  arma::mat M;
  std::vector<std::string> names;
  arma::vec interTemp = arma::zeros(treatbs.n_rows);
  
  for (unsigned int l = 0; l < indexCurrs.size(); ++l) {
    arma::vec indexCurr = indexCurrs[l];
    int i = indexCurr(0);
    int j = indexCurr(1);
    int k = indexCurr(2);

    interTemp = treatSubsamp.col(i) % XSubsamp.col(j) % XSubsamp.col(k);
    
    int q_j = lower_bound(colSizes.begin(), colSizes.end(), j) - colSizes.begin();
    int r_j = 0;
    if (colSizes[q_j] != j) {
      --q_j;
      r_j = j - colSizes[q_j];
    }
    
    int q_k = lower_bound(colSizes.begin(), colSizes.end(), k) - colSizes.begin();
    int r_k = 0;
    if (colSizes[q_j] != k) {
      --q_k;
      r_k = k - colSizes[q_k];
    }
    
    std::string name = treatName + "_bs" + to_string(i) + "_x_" + Xname[q_j] + "_bs" + to_string(r_j) + "_x_" + Xname[q_k] + "_bs" + to_string(r_k);
    names.push_back(name);
    M = arma::join_rows(M, interTemp);
  }
  
  return List::create(Named("cors") = indexCurrs, _["M"] = M, _["names"] = names); //returns a highest correlations, matrix M, variable names
} 

//[[Rcpp::export]]
arma::vec checkcor(arma::mat cors, double thresh) {
  arma::vec v = ones(cors.n_cols); //include all vars initially. This is a bitmask where 1 means to include var and 0 means not to
  for (unsigned int i = 0; i < cors.n_rows; ++i) {
    if (stddev(cors.col(i)) != 0) {
      for (unsigned int j = i; j < cors.n_cols; ++j) {
        if (abs(cors(i, j)) > thresh) {
          v(j) = 0;
        }
      }
    }
    else {
      v(i) = 0;
    }
  }
  return v; //vars marked zero are ones to not include
}

/*//[[Rcpp::export]]
arma::mat gramschmidt(arma::vec y, arma::mat X) { //not finished yet. Right now it just finds the column with the max correlation
 double m = 0;
 double index = -1;
 for (unsigned int i = 0; i < X.n_cols; ++i) {
 arma::vec col = X.col(i);
 if (as_scalar(cor(col, y)) > m) {
    index = i;
    m = as_scalar(cor(col, y));
 }
 }
}*/

//[[Rcpp::export]]
int main() {
  
  //tests
   /*arma::vec a{3,1,4,1,5,9,2,6}; 
   arma::mat b = bs2(a, 3);
   //b.print();
   cout << "\n" << endl;
   arma::mat c = dbs2(a, 3);
   //c.print();
   
   cout << "\n" << endl;
   
   arma::mat d = bsme(a);
   //d.print();
   
   cout << "\n" << endl;
   
   arma::mat e = dbsme(a);
   //e.print();
   
   cout << "\n" << endl;
   
   arma::mat f = makebs(a);
   
   //f.print(); 
   
   arma::mat g{{1,2,3},{4,5,6},{7,8,9}};
   
   cout << "\n" << endl;
   List L = splineBases(g, 3);
   arma::mat h = L["Xbs"];
   h.print();*/
  
  arma::mat X = randn(1000, 5);
  arma::vec treat = randn(1000);
  arma::vec y = 3*randu(1000) + 7;
  List L = splineBasesAndCorrs(1000, 5, X, {"education", "race", "ethnicity", "income", "other"}, y, treat, "treatname", 20); 
  vector<vec> cors = L["cors"];
  for (vec v : cors) {
    cout << v(3) << endl;
  }
  return 0;
}

/*** R
main()
*/
