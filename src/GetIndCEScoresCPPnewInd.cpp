#include <RcppEigen.h>
#include <algorithm>    // to get std::lower_bound
#include <iterator>     // to get std::iterator
#include <limits>       // to get NaN

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
 
Rcpp::List GetIndCEScoresCPPnewInd( const Eigen::Map<Eigen::VectorXd> & yVec, const Eigen::Map<Eigen::VectorXd> & muVec, const Eigen::Map<Eigen::VectorXd> & lamVec, const Eigen::Map<Eigen::MatrixXd> & phiMat, const Eigen::Map<Eigen::MatrixXd> & SigmaYi, const Eigen::Map<Eigen::MatrixXd> & newPhi, const Eigen::Map<Eigen::VectorXd> & newMu ){ 

  // Setting up initial values

  const unsigned int lenyVec = yVec.size();
  const unsigned int lenlamVec = lamVec.size(); 
 
  Eigen::MatrixXd xiVar = Eigen::MatrixXd::Constant(lenlamVec,lenlamVec,std::numeric_limits<double>::quiet_NaN());
  Eigen::MatrixXd xiEst = Eigen::MatrixXd::Constant(lenlamVec,1,std::numeric_limits<double>::quiet_NaN());
  Eigen::MatrixXd fittedY = Eigen::MatrixXd::Constant(lenlamVec,1,std::numeric_limits<double>::quiet_NaN());

  Eigen::MatrixXd LamPhi = lamVec.asDiagonal() * phiMat.transpose();
  Eigen::MatrixXd LamPhiSig = LamPhi * SigmaYi.inverse();

  xiEst = LamPhiSig * (yVec - muVec);
  xiVar = -LamPhi * LamPhiSig.transpose(); 
  xiVar.diagonal() += lamVec;
  fittedY = newMu + newPhi * xiEst;

  return Rcpp::List::create(Rcpp::Named("xiEst") = xiEst,
                            Rcpp::Named("xiVar") = xiVar,
                            Rcpp::Named("fittedY") = fittedY);
}
