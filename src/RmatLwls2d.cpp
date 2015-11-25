#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect and sort

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]


Eigen::MatrixXd RmatLwls2d(
        const Eigen::Map<Eigen::VectorXd> & bw, 
        const std::string kernel_type, 
        const Eigen::Map<Eigen::MatrixXd> & yMat, 
        const Eigen::Map<Eigen::MatrixXd> & wMat, 
        const Eigen::Map<Eigen::VectorXd> & x1In, // Input grid
        const Eigen::Map<Eigen::VectorXd> & x2In, // Input grid
        const Eigen::Map<Eigen::VectorXd> & x1Out,  // Output grid
        const Eigen::Map<Eigen::VectorXd> & x2Out,  // Output grid
        const bool & bwCheck, 
        const bool & transp = true){ 

  // x1Out: out1 (in MATLAB code)
  // x2Out: out2 (in MATLAB code)
  // bwCheck : boolean/ cause the function to simply run the bandwidth check.

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
    //Rcpp::Rcout << "Kernel_type argument was not set correctly; Epanechnikov kernel used." << std::endl;
    Rcpp::warning("Kernel_type argument was not set correctly; Epanechnikov kernel used.");
    KernelName = possibleKernels.find( "epan" )->second;;
  }

  // Start the actual smoother here  
  unsigned int x1OutSize = x1Out.size();
  unsigned int x2OutSize = x2Out.size();
  
  Eigen::MatrixXd mu(x1OutSize, x2OutSize);

  if ( KernelName == 3) { 
    Rcpp::stop("`RmatLwls2d` is slower than `Rmullwlsk` using Gaussian kernel");
  }  

  const double* x1InDat = x1In.data();
  const double* x2InDat = x2In.data();
  const unsigned int x1InSize = x1In.size();
  const unsigned int x2InSize = x2In.size();
  // std::ptrdiff_t* x1UpInd = new std::ptrdiff_t[x1InSize];
  // std::ptrdiff_t* x1LoInd = new std::ptrdiff_t[x1InSize];
  // std::ptrdiff_t* x2UpInd = new std::ptrdiff_t[x2InSize];
  // std::ptrdiff_t* x2LoInd = new std::ptrdiff_t[x2InSize];
  Eigen::VectorXi x1UpInd(x1InSize);
  Eigen::VectorXi x1LoInd(x1InSize);
  Eigen::VectorXi x2UpInd(x2InSize);
  Eigen::VectorXi x2LoInd(x2InSize);

  // Find the x-range (upper and lower ranges are inclusive)
  for (unsigned int i = 0; i < x1OutSize; ++i) {
    const double* uploc = std::lower_bound(x1InDat, x1InDat + x1InSize,
                                     x1Out(i) + bw(0) + 1e-10);
    const double* loloc = std::lower_bound(x1InDat, x1InDat + x1InSize,
                                     x1Out(i) - bw(0) - 1e-10);
    x1UpInd(i) = uploc - x1InDat - 1;
    x1LoInd(i) = loloc - x1InDat;
    if (x1UpInd(i) < x1LoInd(i)) {
      Rcpp::stop("Not enough points in the local window");
    }
  }

  // Find the y-range (upper and lower ranges are inclusive)
  for (unsigned int j = 0; j < x2OutSize; ++j) {
    const double* uploc = std::lower_bound(x2InDat, x2InDat + x2InSize,
                                     x2Out(j) + bw(1) + 1e-10);
    const double* loloc = std::lower_bound(x2InDat, x2InDat + x2InSize,
                                     x2Out(j) - bw(1) - 1e-10);
    x2UpInd(j) = uploc - x2InDat - 1;
    x2LoInd(j) = loloc - x2InDat;
    if (x2UpInd(j) < x2LoInd(j)) {
      Rcpp::stop("Not enough points in the local window");
    }
  }

  // for (int i = 0; i < x1OutSize; ++i) {
    // std::cout << x1LoInd(i) << x1UpInd(i) << std::endl;
  // }
  // for (int j = 0; j < x2OutSize; ++j) {
    // std::cout << x2LoInd(j) << x2UpInd(j) << std::endl;
  // }

  // return(yMat);
  for (unsigned int j = 0; j < x2OutSize; ++j){  
    for (unsigned int i = 0; i < x1OutSize ; ++i){ 

      unsigned int subSize1 = x1UpInd(i) - x1LoInd(i) + 1,
                   subSize2 = x2UpInd(j) - x2LoInd(j) + 1,
                   maxWinSize = subSize1 * subSize2,
                   nnz = 0; // number of non-zero weight points in the
                            // window.
      Eigen::Matrix2Xd lx(2, maxWinSize);
      Eigen::VectorXd ly(maxWinSize);
      Eigen::VectorXd lw(maxWinSize);  
      Eigen::MatrixXd subyMat = yMat.block(x1LoInd(i), x2LoInd(j), 
                                           subSize1, subSize2),
                      subwMat = wMat.block(x1LoInd(i), x2LoInd(j), 
                                           subSize1, subSize2);
      for (unsigned int jw = 0; jw < subSize2; ++jw) {
        for (unsigned int iw = 0; iw < subSize1; ++iw) {
          if (subwMat(iw, jw) > 1e-10) {
            lx(0, nnz) = x1In(iw + x1LoInd(i));
            lx(1, nnz) = x2In(jw + x2LoInd(j));
            ly(nnz) = subyMat(iw, jw);
            lw(nnz) = subwMat(iw, jw);
            ++nnz;
          }
        }
      }
      lx = lx.leftCols(nnz).eval();
      ly = ly.head(nnz).eval();
      lw = lw.head(nnz).eval();


      // check enough points are in the local window 
      unsigned int meter=1;  
      for (unsigned int u =0; u< nnz; ++u) { 
        for (unsigned int t = u + 1; t < nnz; ++t) {
          if ( (lx(0,u) !=  lx(0,t) ) || (lx(1,u) != lx(1,t) ) ) {
            meter++;
          }
        }
        if (meter >= 3) { 
          break; 
        }
      }
   
      //computing weight matrix 
      if (meter >=  3 && !bwCheck) { 
        Eigen::VectorXd temp(nnz);
        Eigen::MatrixXd llx(2, nnz );  
        llx.row(0) = (lx.row(0).array() - x1Out(i))/bw(0);  
        llx.row(1) = (lx.row(1).array() - x2Out(j))/bw(1); 

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
        Eigen::MatrixXd X(nnz ,3);
        X.col(0).setOnes();    
        X.col(1) = lx.row(0).array() - x1Out(i);
        X.col(2) = lx.row(1).array() - x2Out(j); 
        Eigen::LLT<Eigen::MatrixXd> llt_XTWX(X.transpose() * temp.asDiagonal() *X);
        // Probably change LLT to LDLT for increased numerical stability.
        // The solver should stop if the value is NaN. See the HOLE example in gcvlwls2dV2.
        Eigen::VectorXd beta = llt_XTWX.solve(X.transpose() * temp.asDiagonal() * ly);
        mu(i,j)=beta(0); 
      } else if(meter < 3){
        // Rcpp::Rcout <<"The meter value is:" << meter << std::endl;  
        if (bwCheck) {
            Eigen::MatrixXd checker(1,1);
            checker(0,0) = 0.;
            return(checker);
        } else {
          Rcpp::stop("No enough points in local window, please increase bandwidth.");
        }
      }
    }
  }


  if (bwCheck){
     Eigen::MatrixXd checker(1,1); 
     checker(0,0) = 1.; 
     return(checker);
  }
    
  return( mu );
}

 
