#include <vector>
#include <iostream>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma; 

//[[Rcpp::export]]
void OLS(mat X, mat Y) {
    vec beta = inv(trans(X) * X) * (trans(X) * Y);
    vec fitted = X*beta;
    vec errors = abs(Y - fitted);
    
    beta.print();
    fitted.print();
    errors.print();
    
}

//[[Rcpp::export]]
int fib(int x) {
    if (x < 2) {
      return x;
    }
    return fib(x - 1) + fib(x - 2);
}
//[[Rcpp::export]]
int main() {
    
    mat X{{1, 2},{3, 4}};
    vec Y{7, 8};
    OLS(X, Y);
    return 0;
  
}

/***R
main()
*/
