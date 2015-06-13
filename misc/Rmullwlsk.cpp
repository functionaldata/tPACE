#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to read in the kernel name
#include <algorithm>    // to get the intersect and sort

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd Rmullwlsk( const Eigen::Map<Eigen::VectorXd> & bw, const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd> & tPairs, const Eigen::Map<Eigen::MatrixXd> & cxxn, const Eigen::Map<Eigen::VectorXd> & win,  const Eigen::Map<Eigen::VectorXd> & xgrid, const Eigen::Map<Eigen::VectorXd> & ygrid){ 

  // tPairs : xin (in MATLAB code)
  // cxxn : yin (in MATLAB code)
  // xgrid: out1 (in MATLAB code)
  // ygrid: out2 (in MATLAB code)

  // Map the kernel name so we can use switches  
  std::map<std::string,int> possibleKernels; 
  possibleKernels["epan"]    = 1;   possibleKernels["rect"]    = 2;
  possibleKernels["gauss"]   = 3;   possibleKernels["gausvar"] = 4; 
  possibleKernels["quar"]    = 5; 
  int KernelName =  possibleKernels.find(kernel_type )->second; //Set kernel choice


  // Check that we do not have zero weights // Should do a try-catch here
  if ( !(win.all()) ){  // 
    Rcpp::Rcout << "Cases with zero-valued windows are not yet implemented" << std::endl;
    return (tPairs);
  } 

  // Start the actual smoother here  
  unsigned int xgridN = xgrid.size();  
  unsigned int ygridN = ygrid.size();  
  
  Eigen::MatrixXd mu(ygrid.size(),xgrid.size());

  for (unsigned int i = 0; i != 1; ++i){  
    for (unsigned int j = i; j != 3 ; ++j){ 

      //locating local window (LOL)
      std::vector <unsigned int> indx; 
      //if the kernel is not Gaussian
      if ( KernelName != 3) { 
        //construct listX as vectors / size is unknown originally
        std::vector <unsigned int> list1, list2; 
        for (unsigned int y = 0; y != tPairs.cols(); y++){ 
          if ( (tPairs(0,y) >= xgrid(j) - (bw(0)-pow(10,-6)) ) && (tPairs(0,y) <= xgrid(j) + (bw(0)+ pow(10,-6))) ) {
    list1.push_back(y);
          }         
          if ( (tPairs(1,y) >= ygrid(j) - (bw(1)-pow(10,-6)) ) && (tPairs(1,y) <= ygrid(j) + (bw(1)+ pow(10,-6))) ) {
            list2.push_back(y);
          }
        }
        //get intersection 
        std::set_intersection(list1.begin(), list1.begin() + list1.size(), list2.begin(), list2.begin() + list2.size(), std::back_inserter(indx));   
      } else{ // just get the whole deal
        for (unsigned int y = 0; y != tPairs.cols(); ++y){
          indx.push_back(y);
        }
      } 

     for (int yy = 0; yy!= indx.size(); ++yy){   Rcpp::Rcout << indx.at(yy) << ", " ;  }
  Rcpp::Rcout <<  std::endl;
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
      if (meter >=  3) { 
        Eigen::VectorXd temp(indxSize);
        Eigen::MatrixXd llx(2, indxSize );  
        llx.row(0) = (lx.row(0).array() - xgrid(j))/bw(0);  
        llx.row(1) = (lx.row(1).array() - ygrid(i))/bw(1); 

        //define the kernel used 
        switch (KernelName){
          case 1:
         	  Rcpp::Rcout <<" temp.size():"<<std::endl  << temp.size() <<std::endl  ;  
            temp=  ((1-llx.row(0).array().pow(2))*(1- llx.row(1).array().pow(2))).array() * ((9./16)*lw).transpose().array()  ; 
            break;  
          case 2 :
            temp=(lw.array())*.25 ; //this might be wrong.
            break;
          case 3 :  
            temp =  ((-.5*(llx.row(1).array().pow(2))).exp())/(sqrt(2.*M_PI))  *   
              ((-.5*(llx.row(0).array().pow(2))).exp())/(sqrt(2.*M_PI))  *
              (lw.array()); 
            break;
          case 4 :
            temp =  ((((-.5*llx.row(0)).array()).exp())/(sqrt(2.*M_PI))).array() * 
              ((((-.5*llx.row(1)).array()).exp())/(sqrt(2.*M_PI))).array() * 
              (lw.array()) * ((1.5- (((.5) * llx.row(1)).array())).array()) *
              ((1.25- (((.25) * llx.row(0)).array())).array());
            break;
          case 5 :  
            temp =  (lw.array()) * ((1-llx.row(0).array()).array().pow(2)) *
               ((1-llx.row(1).array()).array().pow(2))*(225./256.);
            break;
         } 

         // make the design matrix
         Eigen::MatrixXd X(indxSize ,3);
         X.setOnes();    
         X.col(1) = lx.row(0).array() - xgrid(j);
         X.col(2) = lx.row(1).array() - ygrid(i); 


         //Rcpp::Rcout <<"X : " << std::endl << X << std::endl;  
         // MatrixXd W1(k,k)  = W;   // Something like that should be used
         // Eigen::LLT<Eigen::MatrixXd> llt_XTWX(X.transpose()* temp.asDiagonal() * X);         // Calculate the LLT decomposition
         Eigen::VectorXd beta = (X.transpose()* temp.asDiagonal() * X).inverse() * (X.transpose() * temp.asDiagonal() * ly);  

  
         Rcpp::Rcout <<"llx : " << std::endl << llx << std::endl;  
         //Rcpp::Rcout <<"temp : " << std::endl << temp << std::endl;  
         //Rcpp::Rcout <<"X'WX : " << std::endl << (X.transpose()* temp.asDiagonal()* X) << std::endl;  
         //Eigen::VectorXd beta = llt_XTWX.solve(X.transpose() * temp.asDiagonal() * ly);
         mu(i,j)=beta(0); 
       } else {
         Rcpp::Rcout <<"No enough points in local window, please increase bandwidth." << std::endl;  
         return (tPairs);
       }
    }
  }

  //  This function implements the functionality of mullwlsk

  //Eigen::VectorXd yf = Yst - Xst * betaold;      
  //Eigen::LLT<Eigen::MatrixXd> llt_VY(VY);         // Calculate the LLT decomposition 
        

  // m1=  mu.triangularView<StrictlyUpper>().transpose();
  // m2=  mu.triangularView<Upper>()  ; 
  // m1 = m1+m2;

  //return ( BSold * Zst.transpose() * llt_VY.solve(yf) );    
  return (  mu );
}

 
