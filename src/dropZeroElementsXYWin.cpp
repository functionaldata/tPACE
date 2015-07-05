#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect and sort

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]


Eigen::MatrixXd dropZeroElementsXYWin( const Eigen::Map<Eigen::VectorXd> & win, const Eigen::Map<Eigen::VectorXd> & xin, const Eigen::Map<Eigen::VectorXd> & yin){

  const unsigned int nXGrid = xin.size();

  // Check that we have equal number of readings
  if( nXGrid != yin.size()){
    Rcpp::stop("The input Y-grid does not have the same number of points as input X-grid.");
  }

  if( nXGrid != win.size()){
    Rcpp::stop("The input weight vector does not have the same number of points as input X-grid.");
  }

  unsigned int nZeroElements = std::count(&win[0], &win[nXGrid], 0.);
  
  // Check that we do not have zero weights // Should do a try-catch here
  if( nZeroElements != 0 ){  //
    Eigen::MatrixXd Q(nXGrid - nZeroElements,3);
    unsigned int q = 0;
    for( unsigned int i = 0; i != nXGrid; ++i){
      if ( win[i] != 0 ) {
          Q(q,0) = xin[i];
          Q(q,1) = yin[i];
          Q(q,2) = win[i];
          ++q;
      }
    }
    return( Q );         
  } else {
    Eigen::MatrixXd Q(nXGrid,3);
    Q.col(0) = xin;
    Q.col(1) = yin;
    Q.col(2) = win;
    return( Q );  
  }
}
