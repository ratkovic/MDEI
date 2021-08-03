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




struct Comp { //this is a comparator, used for the heap (priority_queue) in the function below
  public:
    bool operator()(arma::vec a, arma::vec b) {
      return a(3) > b(3);
    }
};

//[[Rcpp::export]]
List splineCorrs(arma::mat XSubsamp, arma::vec ySubsamp, arma::mat treatSubsamp, arma::mat XConstruct, arma::mat treatConstruct,  long long unsigned int a) {
  //a is number of top results, i.e. top 100 or top 300

  priority_queue<arma::vec, std::vector<arma::vec>, Comp> pq;
  arma::vec indexCurr=arma::zeros(4);
  double cor_temp = 0;
  
  for (double i = 0; i < treatSubsamp.n_cols; ++i) {
    for (double j = 0; j < XSubsamp.n_cols - 1; ++j) {
      for (double k = j + 1; k < XSubsamp.n_cols; ++k) {
        arma::vec inter_temp = treatSubsamp.col(i) % XSubsamp.col(j) % XSubsamp.col(k); 
        
        cor_temp = 0;
        if (arma::stddev(inter_temp)>0 ) {
          cor_temp = abs(as_scalar(arma::cor(ySubsamp, inter_temp)));
        }
        //if(!arma::is_finite(cor_temp)) cor_temp =0;
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
// std::vector<std::string> names;
  //arma::vec interTemp = arma::zeros(treatbs.n_rows);
  
  for (unsigned int l = 0; l < indexCurrs.size(); ++l) {
    arma::vec indexCurr = indexCurrs[l];
    int i = indexCurr(0);
    int j = indexCurr(1);
    int k = indexCurr(2);
    
    arma::vec interConstruct = treatConstruct.col(i) % XConstruct.col(j) % XConstruct.col(k);
    // interConstruct = (interConstruct - mean(interConstruct))/stddev(interConstruct);
    
    /*
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
     */
    M = arma::join_rows(M, interConstruct);
  }
  
  return List::create(Named("cors") = indexCurrs, _["M"] = M); //returns a highest correlations, matrix M, variable names
} 


//[[Rcpp::export]]
List namesAndCorrs(arma::mat XSubsamp, arma::vec ySubsamp, arma::mat treatSubsamp, arma::mat XConstruct, arma::mat treatConstruct, arma::mat XConstructDerivative, arma::mat treatConstructDerivative,  long long unsigned int a) {
  //a is number of top results, i.e. top 100 or top 300
  
  priority_queue<arma::vec, std::vector<arma::vec>, Comp> pq;
  arma::vec indexCurr=arma::zeros(4);
  double cor_temp = 0;
  
  for (double i = 0; i < treatSubsamp.n_cols; ++i) {
    for (double j = 0; j < XSubsamp.n_cols - 1; ++j) {
      for (double k = j + 1; k < XSubsamp.n_cols; ++k) {
        arma::vec inter_temp = treatSubsamp.col(i) % XSubsamp.col(j) % XSubsamp.col(k); 
        
        
          cor_temp = abs(as_scalar(arma::cor(ySubsamp, inter_temp)));
          if(!arma::is_finite(cor_temp)) cor_temp = 0;
          
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
  
  arma::mat Msubsamp;
  arma::mat MConstruct;
  arma::mat MConstructDerivative;
  
  for (unsigned int l = 0; l < indexCurrs.size(); ++l) {
    arma::vec indexCurr = indexCurrs[l];
    int i = indexCurr(0);
    int j = indexCurr(1);
    int k = indexCurr(2);
    arma::vec interSubsamp = treatSubsamp.col(i) % XSubsamp.col(j) % XSubsamp.col(k);
    interSubsamp = (interSubsamp - mean(interSubsamp))/stddev(interSubsamp);
    
    arma::vec interConstruct = treatConstruct.col(i) % XConstruct.col(j) % XConstruct.col(k);
    interConstruct = (interConstruct - mean(interConstruct))/stddev(interConstruct);
    
    arma::vec interConstructDerivative = treatConstructDerivative.col(i) % XConstructDerivative.col(j) % XConstructDerivative.col(k);
    interConstructDerivative = (interConstructDerivative - mean(interConstructDerivative))/stddev(interConstructDerivative);
    
    Msubsamp = arma::join_rows(Msubsamp, interSubsamp);
    MConstruct = arma::join_rows(MConstruct, interConstruct);
    MConstructDerivative = arma::join_rows(MConstructDerivative, interConstructDerivative);
  }
  
  return List::create(Named("cors") = indexCurrs, _["Msubamp"] = Msubsamp, _["MConstruct"] = MConstruct, _["MConstructDerivative"] = MConstructDerivative); //returns a highest correlations, matrix M, variable names
} 