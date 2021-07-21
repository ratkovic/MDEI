// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <queue>
#include <vector>
#include <time.h>
#include <map>

using namespace arma;
using namespace Rcpp;
using namespace std;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' Check Spearman correlations between interactions in X and treatment
//' 
//' @param y A vector of outcomes.
//' @param X A matrix of spline bases.
//' @param alpha.schedule The prior on lambda
//' @export

//[[Rcpp::export]]
arma::vec myrank(arma::vec x) {
  arma::vec sorted = arma::sort(x);
  arma::vec rank = arma::zeros(x.size());
  unordered_map<int, int> umap;
  unordered_map<int, int> counts;
  
  for (unsigned int i = 0; i < sorted.size(); ++i) {
    umap[sorted[i]] += i + 1;
    counts[sorted[i]]++;
  }
  
  for (unsigned int i = 0; i < sorted.size(); ++i) {
    rank(i) = double(umap[x(i)]/(counts[x(i)] + 0.0));
  }
  return rank;
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
List namesAndCorrs(arma::mat XSubsamp, std::vector<std::string> Xnames, arma::vec ySubsamp, std::vector<int> colSizes, arma::mat treatSubsamp, arma::mat XConstruct, arma::mat treatConstruct, std::vector<std::string> treatNames, long long unsigned int a) {
  //a is number of top results, i.e. top 100 or top 300

  priority_queue<arma::vec, std::vector<arma::vec>, Comp> pq;
  arma::vec indexCurr=arma::zeros(4);

  for (double i = 0; i < treatSubsamp.n_cols; ++i) {
    for (double j = 0; j < XSubsamp.n_cols - 1; ++j) {
      for (double k = j + 1; k < XSubsamp.n_cols; ++k) {
        arma::vec inter_temp = treatSubsamp.col(i) % XSubsamp.col(j) % XSubsamp.col(k); 
        
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
  
  for (unsigned int l = 0; l < indexCurrs.size(); ++l) {
    arma::vec indexCurr = indexCurrs[l];
    int i = indexCurr(0);
    int j = indexCurr(1);
    int k = indexCurr(2);
    
    arma::vec interConstruct = treatConstruct.col(i) % XConstruct.col(j) % XConstruct.col(k);
    interConstruct = (interConstruct - mean(interConstruct))/stddev(interConstruct);
    
    std::string name = treatNames[i] + "_bs" + to_string(i) + "_x_" + Xnames[j] + "_bs" + to_string(j) + "_x_" + Xnames[k] + "_bs" + to_string(k);
    names.push_back(name);
    M = arma::join_rows(M, interConstruct);
  }
  
  return List::create(Named("cors") = indexCurrs, _["M"] = M, _["names"] = names); //returns a highest correlations, matrix M, variable names
} 
/*
//[[Rcpp::export]]
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
/*
//[[Rcpp::export]]
int main() {
  
  //tests
   arma::vec a{3,1,4,1,5,9,2,6,5,3};

   arma::mat b = bs2(a, 3);
   //b.print();
   cout << "\n" << endl;
   arma::mat c = dbs2(a, 3);
   //c.print();
   
   cout << "\n" << endl;
   
   arma::mat d = bsme(a);
   d.print();
   
   cout << "\n" << endl 
   
   arma::mat e = dbsme(a);
   //e.print();
   
   cout << "\n" << endl;
   
   arma::mat f = makebs(a);
   
   //f.print(); 
   
   arma::mat g{{1,2,3},{4,5,6},{7,8,9}};
   
   cout << "\n" << endl;
   List L = splineBases(g, 3);
   arma::mat h = L["Xbs"];
   h.print(); 
  
  arma::mat X = randn(10, 5);
  arma::vec treat = randn(10);
  arma::vec y = 3*randu(10) + 7;
  arma::vec replaceme = {1, 1, 1, 1, 2, 1, 2, 2, 2, 2};
  arma::vec alphas{10,9,8,7,6,5,4,3,2,1};
  List L = namesAndCorrs(X, {"education", "race", "ethnicity", "income", "other"}, y, replaceme, treat, "treatname", 20); 
  arma::mat M = L["M"];
  M.print();
  List L2 = GCV(y, M, alphas, 1e-4);
  double gcv = L2["GCV"];
  cout << gcv << endl;
  return 0;
}*/