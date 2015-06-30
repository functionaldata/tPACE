#include <RcppEigen.h>
#include <algorithm>    // to get std::lower_bound
#include <iterator>     // to get std::iterator

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::VectorXd interp2lin( const Eigen::Map<Eigen::VectorXd> & xin, const Eigen::Map<Eigen::VectorXd> & yin, const Eigen::Map<Eigen::VectorXd> & zin, const Eigen::Map<Eigen::VectorXd> & xou, const Eigen::Map<Eigen::VectorXd> & you){ 

  // Setting up initial values

  const unsigned int nXGrid = xin.size();
  const unsigned int nYGrid = yin.size();
  const unsigned int nKnownPoints = nXGrid * nYGrid;
  const unsigned int nUnknownPoints = xou.size();
 
  Eigen::VectorXd result(nUnknownPoints);

  // Preliminary checks

  if ( nXGrid != nYGrid ){
    Rcpp::stop("Input Y-grid does not have the same number of points as Input X-grid.");
  } else if ( nKnownPoints != zin.size() ) {  
    Rcpp::stop("Input Z-grid does not have the same number of points as the product of #Input Y-grid times #Input X-grid.");
  } else if ( nUnknownPoints != you.size() ){
    Rcpp::stop("Output Y-grid does not have the same number of points as Output X-grid.");
  } else if ( xin.minCoeff() >  xou.minCoeff() ){
    Rcpp::stop("Output X-grid  is outside the lower range of the input X-grid.");
  } else if ( yin.minCoeff() >  you.minCoeff() ){
    Rcpp::stop("Output X-grid  is outside the lower range of the input X-grid.");
  } else if ( xin.maxCoeff() <  xou.maxCoeff() ){
    Rcpp::stop("Output X-grid  is outside the upper ragne of the input X-grid.");
  } else if ( yin.maxCoeff() <  you.maxCoeff() ){
    Rcpp::stop("Output X-grid  is outside the upper ragne of the input X-grid.");
  } 

  // The actual interpolation routine
  // Check: https://en.wikipedia.org/wiki/Bilinear_interpolation#Alternative_algorithm
  // for a quick reference

  Eigen::RowVector4d fq;
  Eigen::Vector2d xa;
  Eigen::Vector2d ya;
  Eigen::Vector4d za;
  Eigen::Matrix4d A;
  // Column 1
  A.setOnes(); 

  for (unsigned int u = 0; u !=nUnknownPoints; ++u){ 
  
    // Find the appropriate x coordinates/save them in xa (2-by-1)
    // Get iterator pointing to the first element which is not less than xou(u)
    const double* x1 = std::lower_bound(&xin[0], &xin[nXGrid], xou(u)); 
    xa(1) = *x1;
    xa(0) = *--x1;
  
    // I feel there is a bug around here. But I cannot smoke it out.
  
    // Find the appropriate y coordinates/save them in ya (2-by-1)
    // Get iterator pointing to the first element which is not less than za(u)
    const double* y1 = std::lower_bound(&yin[0], &yin[nYGrid], you(u)); 
    ya(1) = *y1;
    ya(0) = *--y1;
  
    za(0) =  zin( (std::find(&xin[0], &xin[nXGrid], xa(0)) -&xin[0]) * nXGrid + (std::find(&yin[0], &yin[nXGrid], ya(0)) -&yin[0]));
    za(1) =  zin( (std::find(&xin[0], &xin[nXGrid], xa(0)) -&xin[0]) * nXGrid + (std::find(&yin[0], &yin[nXGrid], ya(1)) -&yin[0]));
    za(2) =  zin( (std::find(&xin[0], &xin[nXGrid], xa(1)) -&xin[0]) * nXGrid + (std::find(&yin[0], &yin[nXGrid], ya(0)) -&yin[0]));
    za(3) =  zin( (std::find(&xin[0], &xin[nXGrid], xa(1)) -&xin[0]) * nXGrid + (std::find(&yin[0], &yin[nXGrid], ya(1)) -&yin[0]));

    // Column 2  
    A(0,1) = xa(0);   A(1,1) = xa(0);
    A(2,1) = xa(1);   A(3,1) = xa(1);
    // Column 3
    A(0,2) = ya(0);   A(1,2) = ya(1);
    A(2,2) = ya(0);   A(3,2) = ya(1);
    // Column 4
    A(0,3) = ya(0) * xa(0);   A(1,3) = ya(1) * xa(0);
    A(2,3) = ya(0) * xa(1);   A(3,3) = ya(1) * xa(1);
  
    fq << 1 , xou(u), you(u), xou(u)*you(u);
  
    // Rcpp::Rcout << "xa: " << xa.transpose() << ", ya: " << ya.transpose()  << ", za: " << za.transpose() << ", fq: " << fq << std::endl;
  
    result(u) = fq * A.colPivHouseholderQr().solve(za);
  
  }
  
  return ( result ); 
}
