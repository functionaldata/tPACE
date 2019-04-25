#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector Rcppsort(NumericVector v) {
  NumericVector sv(clone(v));
  std::sort(sv.begin(), sv.end());
  return sv;
}

 
