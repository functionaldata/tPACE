#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect and sort

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::VectorXd Rrotatedmullwlsk( const Eigen::Map<Eigen::VectorXd> & bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd> & tPairs, const Eigen::Map<Eigen::MatrixXd> & cxxn, const Eigen::Map<Eigen::VectorXd> & win,  const Eigen::Map<Eigen::MatrixXd> & xygrid, const unsigned int npoly, const bool & bwCheck){ 

  // tPairs : xin (in MATLAB code)
  // cxxn : yin (in MATLAB code)
  // xygrid: d (in MATLAB code)
  // npoly: redundant?

  const double invSqrt2pi=  1./(sqrt(2.*M_PI));

  // Map the kernel name so we can use switches  
  std::map<std::string,int> possibleKernels; 
  possibleKernels["epan"]    = 1;   possibleKernels["rect"]    = 2;
  possibleKernels["gauss"]   = 3;   possibleKernels["gausvar"] = 4; 
  possibleKernels["quar"]    = 5; 
   
  // The following test is here for completeness, we mightwant to move it up a 
  // level (in the wrapper) in the future. 

  // If the kernel_type key exists set KernelName appropriately
  int KernelName = 0;
  if ( possibleKernels.count( kernel_type ) != 0){ 
    KernelName = possibleKernels.find( kernel_type )->second; //Set kernel choice
  } else {
  // otherwise use "epan"as the kernel_type 
  // Rcpp::Rcout << "Kernel_type argument was not set correctly; Epanechnikov kernel used." << std::endl;
    Rcpp::warning("Kernel_type argument was not set correctly; Epanechnikov kernel used.");
    KernelName = possibleKernels.find( "epan" )->second;;
  }

  // Check that we do not have zero weights // Should do a try-catch here
  // Again this might be best moved a level-up. 
  if ( !(win.all()) ){  // 
    Rcpp::Rcout << "Cases with zero-valued windows are not yet implemented" << std::endl;
    return (tPairs);
  } 

  Eigen::Matrix2d RC;  // Rotation Coordinates
  RC << 1, -1, 1, 1; 
  RC *= sqrt(2.)/2.; 

  Eigen::MatrixXd rtPairs = RC * tPairs;
  Eigen::MatrixXd rxygrid = RC * xygrid;
  
  unsigned int xygridN = rxygrid.cols();  
  
  Eigen::VectorXd mu(xygridN);
  mu.setZero();     

  for (unsigned int i = 0; i !=  xygridN; ++i){   

    //locating local window (LOL) (bad joke)
    std::vector <unsigned int> indx; 
    //if the kernel is not Gaussian or Gaussian-like
    if ( KernelName != 3 && KernelName != 4 ) { 
      //construct listX as vectors / size is unknown originally
      std::vector <unsigned int> list1, list2; 
      for (unsigned int y = 0; y != tPairs.cols(); y++){ 
        if ( std::abs( rtPairs(0,y) - rxygrid(0,i) ) <= bw(0)  ) {
          list1.push_back(y);
        }         
        if ( std::abs( rtPairs(1,y) - rxygrid(1,i) ) <= bw(1)  ) {
          list2.push_back(y);
        }
      }
      
    //get intersection between the two lists 
    std::set_intersection(list1.begin(), list1.begin() + list1.size(), list2.begin(), list2.begin() + list2.size(), std::back_inserter(indx));   
  
    } else { // just get the whole deal
      for (unsigned int y = 0; y != tPairs.cols(); ++y){
        indx.push_back(y);
      }
    }   

    unsigned int indxSize = indx.size();
    Eigen::VectorXd lw(indxSize);  
    Eigen::VectorXd ly(indxSize);
    Eigen::MatrixXd lx(2,indxSize);
    for (unsigned int u = 0; u !=indxSize; ++u){ 
      lx.col(u) = rtPairs.col(indx[u]); 
      lw(u) = win(indx[u]); 
      ly(u) = cxxn(indx[u]); 
    }


    if (ly.size()>=npoly+1 && !bwCheck ){

      //computing weight matrix 
      Eigen::VectorXd temp(indxSize);
      Eigen::MatrixXd llx(2, indxSize );  
      llx.row(0) = (lx.row(0).array() - rxygrid(0,i))/bw(0);  
      llx.row(1) = (lx.row(1).array() - rxygrid(1,i))/bw(1); 
 
      //define the kernel used 
      switch (KernelName){
        case 1: // Epan
          temp=  ((1-llx.row(0).array().pow(2))*(1- llx.row(1).array().pow(2))).array() * 
                 ((9./16)*lw).transpose().array(); 
          break;  
        case 2 : // Rect
          temp=(lw.array())*.25 ; 
          break;
      case 3 : // Gauss
          temp = ((-.5*(llx.row(1).array().pow(2))).exp()) * invSqrt2pi  *   
                 ((-.5*(llx.row(0).array().pow(2))).exp()) * invSqrt2pi  *
                 (lw.transpose().array()); 
          break;
      case 4 : // GausVar
          temp = (lw.transpose().array()) * 
                 ((-0.5 * llx.row(0).array().pow(2)).array().exp() * invSqrt2pi).array() *
                 ((-0.5 * llx.row(1).array().pow(2)).array().exp() * invSqrt2pi).array() * 
                  (1.25 - (0.25 * (llx.row(0).array().pow(2))).array())  * 
                  (1.50 - (0.50 * (llx.row(1).array().pow(2))).array()); 
          break;
        case 5 :  // Quar
          temp = (lw.transpose().array()) * 
                 ((1.-llx.row(0).array().pow(2)).array().pow(2)).array() *
                 ((1.-llx.row(1).array().pow(2)).array().pow(2)).array() * (225./256.);
          break;
      } 
      
      // make the design matrix
      Eigen::MatrixXd X(indxSize ,3);
      X.setOnes();    
      X.col(1) = (lx.row(0).array() - rxygrid(0,i)).array().pow(2);
      X.col(2) = (lx.row(1).array() - rxygrid(1,i)); 
      Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(X.transpose() * temp.asDiagonal() *X);
      Eigen::VectorXd beta = ldlt_XTWX.solve(X.transpose() * temp.asDiagonal() * ly);
      mu(i)=beta(0); 
    } else if ( ly.size() == 1 && !bwCheck) { // Why only one but not two is handled?
      mu(i) = ly(0);
    } else if ( ly.size() != 1 && (ly.size() < npoly+1) ) {
      if ( bwCheck ){
        Eigen::VectorXd checker(1); 
        checker(0) = 0.; 
        return(checker);
      } else {  
        Rcpp::stop("No enough points in local window, please increase bandwidth.");
      }
    }
  } 
  
  if (bwCheck){
     Eigen::VectorXd checker(1); 
     checker(0) = 1.; 
     return(checker);
  }
  
  return ( mu ); 
}

 
