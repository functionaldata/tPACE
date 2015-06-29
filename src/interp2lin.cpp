#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect and sort
#include <iterator>     // std::iterator, std::input_iterator_tag

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::VectorXd interp2lin( const Eigen::Map<Eigen::VectorXd> & xin, const Eigen::Map<Eigen::VectorXd> & yin, const Eigen::Map<Eigen::VectorXd> & zin, const Eigen::Map<Eigen::VectorXd> & xou, const Eigen::Map<Eigen::VectorXd> & you){ 

  // Setting up initial values

  const unsigned int nKnownPoints = xin.size();
  const unsigned int nUnknownPoints = xou.size();
 
  Eigen::VectorXd result(nUnknownPoints);

  // Preliminary checks

  if ( nKnownPoints != yin.size() ){
    Rcpp::stop("Input Y-grid does not have the same number of points as Input X-grid.");
  } else if ( nKnownPoints != zin.size() ) {  
    Rcpp::stop("Input Z-grid does not have the same number of points as Input X-grid.");
  } else if ( nUnknownPoints != you.size() ){
    Rcpp::stop("Output Y-grid does not have the same number of points as Output X-grid.");
  } else if ( xin.minCoeff() >  xou.minCoeff() ){
    Rcpp::stop("Output X-grid  is outside the lower ragne of the input X-grid.");
  } else if ( yin.minCoeff() >  you.minCoeff() ){
    Rcpp::stop("Output X-grid  is outside the lower ragne of the input X-grid.");
  } else if ( xin.maxCoeff() <  xou.maxCoeff() ){
    Rcpp::stop("Output X-grid  is outside the upper ragne of the input X-grid.");
  } else if ( yin.maxCoeff() <  you.maxCoeff() ){
    Rcpp::stop("Output X-grid  is outside the upper ragne of the input X-grid.");
  } 

  // The actual interpolation routine

  // Make sorted copies of xin & xou

  Eigen::VectorXd xinput = xin;
  std::sort( &xinput[0], &xinput[nKnownPoints]);

  Eigen::VectorXd yinput = yin;
  std::sort( &yinput[0], &yinput[nKnownPoints]);

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
  double* x1 = std::lower_bound(&xinput[0], &xinput[nKnownPoints], xou(u)); 
  xa(1) = *x1;
  xa(0) = *--x1;

// I suspect there is a bug around here.

  // Find the appropriate y coordinates/save them in ya (2-by-1)
  // Get iterator pointing to the first element which is not less than za(u)
  double* y1 = std::lower_bound(&yinput[0], &yinput[nKnownPoints], you(u)); 
  ya(1) = *y1;
  ya(0) = *--y1;

  unsigned int foundAll = 0;

  // Find the associated z values/save them in za (4-by-1)
  for (unsigned int i = 0; i !=nKnownPoints; ++i){ 
    if (xa(0) == xin(i) && ya(0) == yin(i)) {
      za(0) = zin(i);
      foundAll++; 
      continue;
    }

    if (xa(0) == xin(i) && ya(1) == yin(i)) {
      za(1) = zin(i); 
      foundAll++;
      continue;
    }

    if (xa(1) == xin(i) && ya(0) == yin(i)) {
      za(2) = zin(i); 
      foundAll++;
      continue;
    }
     
    if (xa(1) == xin(i) && ya(1) == yin(i)) {
      za(3) = zin(i);  
      foundAll++; 
      continue;
    }

    if ( foundAll == 4 ) {
      break;
    } 
  } 

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
