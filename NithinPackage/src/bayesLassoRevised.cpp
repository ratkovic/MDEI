// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace arma;
using namespace std;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

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
Rcpp::List bayesLasso(vec y, mat X, vec alpha, double tol) { //tol is 1e-6
  
  int n = X.n_rows;
  int p = X.n_cols;
  vec Etausqinv = ones(p);
  vec Ewtsqtausq = ones(p);
  mat XpX = X.t()*X;
  mat Xpy = X.t()*y;
  mat XpXsolve = XpX;
  vec fits = zeros(n);
  double lambda_0 = 1;
  
  double sdy = stddev(y);
  double conv = 1;
  double edf = 0;
  double GCV = 0;
  
  vec beta = zeros(p);
  vec beta_last = beta;
  
  int iters = 500; //I like to put constants like this as their own variables, so we can change the value later easily
  
  for (int i = 0; i < iters; ++i) {
    
    if (conv/sdy > tol) {
      for (int j = 1; j < p; ++j) {
        XpXsolve(j, j) = XpX(j, j) + Etausqinv(j) + 1e-6;    
      }
      beta_last = beta;
      uvec update_ind = find(abs(beta) > tol*sdy);
      beta(update_ind) = solve(XpXsolve.submat(update_ind,update_ind),Xpy.rows(update_ind));
      fits = X*beta;
      for (int j= 1; j < p; ++j) {
        Ewtsqtausq(j)  = abs(beta(j))/(lambda_0*sdy)+pow(lambda_0,-2);
        Etausqinv(j)  = lambda_0/abs(beta(j))*sdy;
      }
      Etausqinv(0) = 0;
      Ewtsqtausq(0) = 0;
      
      lambda_0 = sqrt((alpha(0) - 1)/(sum(Ewtsqtausq)/2 + 1));
      //sigma_sq = (sum((y - fits) % (y - fits)) + sum(beta % beta % Etausqinv/2))/(n/2 + p/2 + 1);
      
      conv = max(abs(beta - beta_last));
    }
    uvec update_ind = find(abs(beta) > tol*sdy );
    edf = trace(XpX.submat(update_ind,update_ind)*pinv(XpXsolve.submat(update_ind,update_ind)));
    double den = (n-log(n)/2*edf);
    GCV = sum((y-fits)%(y-fits))/(den*den);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("intercept") = beta(0),
    Rcpp::Named("coefficients") = beta.rows(1,p-1),
    Rcpp::Named("fitted.values") = fits,
    Rcpp::Named("GCV") = GCV,
    Rcpp::Named("Etausqinv") = Etausqinv
  );
}
//[[Rcpp::export]]
int main() {
    vec y{3, 1, 4, 1};
    mat X{{2, 7, 1},{8, 2, 8}, {1, 8, 2}, {8, 4, 5}};
    
    vec alpha{2};
    double tol = 1e-6;
    
    print(bayesLasso(y, X, alpha, tol));
    
    return 0;
}

/*** R
main()
*/
