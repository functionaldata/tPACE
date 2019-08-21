#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double RCPPvar(const Rcpp::NumericVector X){
 return ( var(X) ) ;
}
