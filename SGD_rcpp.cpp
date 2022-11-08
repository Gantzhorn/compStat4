#include <Rcpp.h>


// [[Rcpp::depends(dqrng]]

// [[Rcpp::export]]
Rcpp::IntegerVector cppsample(int n){
  Rcpp::IntegerVector I = Rcpp::sample(n, n);
  return I;
}
