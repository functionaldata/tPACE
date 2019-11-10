#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect, sort, lower_bound, upper_bound
// #include <gperftools/profiler.h>

typedef std::pair<double, unsigned int> valIndPair;
bool comparePair(const valIndPair& l, const valIndPair& r) {
  return(l.first < r.first);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd RmullwlskUniversal( const Eigen::Map<Eigen::VectorXd> & bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd> & tPairs, const Eigen::Map<Eigen::MatrixXd> & cxxn, const Eigen::Map<Eigen::VectorXd> & win,  const Eigen::Map<Eigen::VectorXd> & xgrid, const Eigen::Map<Eigen::VectorXd> & ygrid, const bool & bwCheck, const bool & autoCov){ 
// Assumes the first row of tPairs is sorted in increasing order.
  // tPairs : xin (in MATLAB code)
  // cxxn : yin (in MATLAB code)
  // xgrid: out1 (in MATLAB code)
  // ygrid: out2 (in MATLAB code)
  // bwCheck : boolean / cause the function to simply run the bandwidth check. //To be depreciated
  // autoCov : boolean / cause the function to return the autocovariance.

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

  // Check that we do not have zero weights // Should do a try-catch here
  // Again this might be best moved a level-up. 
  if ( !(win.all()) ){  // 
    Rcpp::Rcout << "Cases with zero-valued windows are not yet implemented" << std::endl;
    return (tPairs);
  } 

  // ProfilerStart("sort.log");
  // Start the actual smoother here  
  const unsigned int xgridN = xgrid.size();  
  const unsigned int ygridN = ygrid.size();  
  const unsigned int n = tPairs.cols();
  
  // For sorted x1
  Eigen::VectorXd x1(tPairs.row(0).transpose());
  const double* tDat = x1.data();
  
  Eigen::MatrixXd mu(xgrid.size(), ygrid.size());
  mu.setZero();    


  for (unsigned int i = 0; i != xgridN; ++i) {  
    const double xl = xgrid(i) - bw(0) - 1e-6, 
                 xu = xgrid(i) + bw(0) + 1e-6;

    unsigned int indl = std::lower_bound(tDat, tDat + n, xl) - tDat,  
                 indu = std::upper_bound(tDat, tDat + n, xu) - tDat;

    // sort the y index
    std::vector<valIndPair> yval(indu - indl);
    for (unsigned int k = 0; k < yval.size(); ++k){
      yval[k] = std::make_pair(tPairs(1, k + indl), k + indl);
    }
    std::sort<std::vector<valIndPair>::iterator>(yval.begin(), yval.end(), comparePair);

    std::vector<valIndPair>::iterator ylIt = yval.begin(), 
                                      yuIt = yval.begin();

    for (unsigned int j = !autoCov? 0: i; j != ygridN; ++j) { 
      const double yl = ygrid(j) - bw(1) - 1e-6, 
                   yu = ygrid(j) + bw(1) + 1e-6;

      //locating local window (LOL) (bad joke)
      std::vector <unsigned int> indx; 
      
      //if the kernel is not Gaussian
      if ( KernelName != 3) { 
      // Search the lower and upper bounds increasingly.
        ylIt = std::lower_bound(ylIt, yval.end(), valIndPair(yl, 0), comparePair);
        yuIt = std::upper_bound(yuIt, yval.end(), valIndPair(yu, 0), comparePair);

        // The following works nice for the Gaussian 
        //  but for very small samples it complains  
        //} else {
        //  ylIt = yval.begin();
        //  yuIt = yval.end();
        //}

        for (std::vector<valIndPair>::iterator y = ylIt; y != yuIt; ++y){ 
          indx.push_back(y->second);
        } 
      } else { //When we finally get c++11 we will use std::iota
        for( unsigned int y = 0; y != n; ++y){
          indx.push_back(y);
        }
      }

      unsigned int indxSize = indx.size();
      Eigen::VectorXd lw(indxSize);  
      Eigen::VectorXd ly(indxSize);
      Eigen::MatrixXd lx(2,indxSize);

      for (unsigned int u = 0; u !=indxSize; ++u){ 
        lx.col(u) = tPairs.col(indx[u]); 
        lw(u) = win(indx[u]); 
        ly(u) = cxxn(indx[u]); 
      }

      // check enough points are in the local window 
      unsigned int meter=1;  
      for (unsigned int u =0; u< indxSize; ++u) { 
        for (unsigned int t = u + 1; t < indxSize; ++t) {
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
        Eigen::VectorXd temp(indxSize);
        Eigen::MatrixXd llx(2, indxSize );  
        llx.row(0) = (lx.row(0).array() - xgrid(i))/bw(0);  
        llx.row(1) = (lx.row(1).array() - ygrid(j))/bw(1); 

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
        X.col(1) = lx.row(0).array() - xgrid(i);
        X.col(2) = lx.row(1).array() - ygrid(j); 
        Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(X.transpose() * temp.asDiagonal() *X);
        // The solver should stop if the value is NaN. See the HOLE example in gcvlwls2dV2.
        Eigen::VectorXd beta = ldlt_XTWX.solve(X.transpose() * temp.asDiagonal() * ly);
        mu(i,j)=beta(0);  
      } 
      else if(meter < 3){
        // // Rcpp::Rcout <<"The meter value is:" << meter << std::endl;  
        // if (bwCheck) {
            // Eigen::MatrixXd checker(1,1);
            // checker(0,0) = 0.;
            // return(checker);
        // } else {
             Rcpp::stop("No enough points in local window, please increase bandwidth using userBwCov.");
        // }
      }
    }
  }

  if (bwCheck){
     Eigen::MatrixXd checker(1,1); 
     checker(0,0) = 1.; 
     return(checker);
  } 
  if(autoCov){
       return ( Eigen::MatrixXd(mu.triangularView<Eigen::StrictlyUpper>().transpose()) + Eigen::MatrixXd(mu.triangularView<Eigen::Upper>()));
  } else {
    // ProfilerStop();
    return ( mu ); 
  }
}

 
