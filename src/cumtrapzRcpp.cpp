
#include <Rcpp.h> 
using namespace Rcpp;

template <class ForwardIterator> bool is_sorted (ForwardIterator first, ForwardIterator last)
{
  if (first==last) return true;
  ForwardIterator next = first;
  while (++next!=last) {
    if (*next<*first)
      return false;
    ++first;
  }
  return true;
}
//' Cumulative Trapezoid Rule Numerical Integration
//' 
//' Cumulative Trapezoid Rule Numerical Integration using Rcpp
//' @param X Sorted vector of X values
//' @param Y Vector of Y values.
//' @export 
// [[Rcpp::export]]
Rcpp::NumericVector cumtrapzRcpp(const Rcpp::NumericVector X,const Rcpp::NumericVector Y){   

  // Basic check
  if( Y.size() != X.size()){
    Rcpp::stop("The input Y-grid does not have the same number of points as input X-grid.");
  }
  if(is_sorted(X.begin(),X.end())){
    Rcpp::NumericVector  ctrapzsum(X.size()); 
    ctrapzsum[0] = 0.0;    
    for (unsigned int ind = 0; ind !=  X.size()-1; ++ind){
      ctrapzsum[ind+1] = 0.5 * (X[ind + 1] - X[ind]) *(Y[ind] + Y[ind + 1]) + ctrapzsum[ind];  
    }
    return ctrapzsum;
  } else {
    Rcpp::stop("The input X-grid is not sorted.");
    return 1;
  }
}
