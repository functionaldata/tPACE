#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::MatrixXd calc_muB( const Eigen::Map<Eigen::MatrixXd> & BSold, const Eigen::Map<Eigen::MatrixXd> & VY, const Eigen::Map<Eigen::MatrixXd> & Xst, const Eigen::Map<Eigen::MatrixXd> & Zst, const Eigen::Map<Eigen::VectorXd> & Yst,  const Eigen::Map<Eigen::VectorXd> & betaold){ 

  //  This function implements:
  //  as.vector(BSold %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]]-X.st[[i]]%*%betaold)))

  Eigen::VectorXd yf = Yst - Xst * betaold;      
  Eigen::LLT<Eigen::MatrixXd> llt_VY(VY);         // Calculate the LLT decomposition 
        
  return ( BSold * Zst.transpose() * llt_VY.solve(yf) );    
 
}
