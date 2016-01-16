#include <Rcpp.h> 
#include <limits>       // to get NaN
using namespace Rcpp;

template <class iter> bool is_sorted (iter begin, iter end)
{
  if (begin==end) return true;
  iter next = begin;
  while (++next!=end) {
    if (*next<*begin)
      return false;
    ++begin;
  }
  return true;
}

// [[Rcpp::export]]
double trapzRcpp(const Rcpp::NumericVector X, const Rcpp::NumericVector Y){   

  if( Y.size() != X.size()){
    Rcpp::stop("The input Y-grid does not have the same number of points as input X-grid.");
  }
  if(is_sorted(X.begin(),X.end())){
    double trapzsum = 0; 
    for (unsigned int ind = 0; ind !=  X.size()-1; ++ind){
      trapzsum += 0.5 * (X[ind + 1] - X[ind]) *(Y[ind] + Y[ind + 1]); 
    }
    return trapzsum;
  } else {
    Rcpp::stop("The input X-grid is not sorted.");
    return  std::numeric_limits<double>::quiet_NaN();
  }
  return  std::numeric_limits<double>::quiet_NaN();
}
