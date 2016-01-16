#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double RCPPmean(const Rcpp::NumericVector X){
 return ( mean(X) ) ;
}
