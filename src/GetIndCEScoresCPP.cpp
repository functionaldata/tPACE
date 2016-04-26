#include <RcppEigen.h>
#include <algorithm>    // to get std::lower_bound
#include <iterator>     // to get std::iterator
#include <limits>       // to get NaN

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
 

Rcpp::List GetIndCEScoresCPP( const Eigen::Map<Eigen::VectorXd> & yVec, const Eigen::Map<Eigen::VectorXd> & muVec, const Eigen::Map<Eigen::VectorXd> & lamVec, const Eigen::Map<Eigen::MatrixXd> & phiMat, const Eigen::Map<Eigen::MatrixXd> & SigmaYi){ 

  // Setting up initial values

  const unsigned int lenyVec = yVec.size();
  const unsigned int lenlamVec = lamVec.size(); 
 
  Eigen::MatrixXd xiVar = Eigen::MatrixXd::Constant(lenlamVec,lenlamVec,std::numeric_limits<double>::quiet_NaN());
  Eigen::MatrixXd xiEst = Eigen::MatrixXd::Constant(lenlamVec,1,std::numeric_limits<double>::quiet_NaN());
  Eigen::MatrixXd fittedY = Eigen::MatrixXd::Constant(lenlamVec,1,std::numeric_limits<double>::quiet_NaN());
   
  Eigen::MatrixXd LamPhi = lamVec.asDiagonal() * phiMat.transpose();  
  Eigen::LDLT<Eigen::MatrixXd> ldlt_SigmaYi(SigmaYi);
  xiEst = LamPhi * ldlt_SigmaYi.solve(yVec - muVec) ;       // LamPhiSig * (yVec - muVec); 
  xiVar = -LamPhi * ldlt_SigmaYi.solve(LamPhi.transpose()); // LamPhiSig.transpose(); 
  
  xiVar.diagonal() += lamVec;
  fittedY = muVec + phiMat * xiEst;

  return Rcpp::List::create(Rcpp::Named("xiEst") = xiEst,
                            Rcpp::Named("xiVar") = xiVar,
                            Rcpp::Named("fittedY") = fittedY);
}
