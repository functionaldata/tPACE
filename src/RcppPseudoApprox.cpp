#include <Rcpp.h>
#include <RcppEigen.h>
 
// [[Rcpp::depends(RcppEigen)]]

float LinearInterpolation ( const Eigen::Map<Eigen::VectorXd> & X ,  const Eigen::Map<Eigen::VectorXd> & Y, float X_PointOfInterest){
  //Produce Y_point_of_interest given X,Y and target X_point_of_interest
  //X : vector containing the X variables of the interpolant 
  //Y : vector containing the Y variables of the interpolant 
  //PointOfInterest : Point of X to estimate the new point of Y
  
  float   xk, xkp1, yk, ykp1 = 0;  //Points adjecent to the point of interpolation
  if ( X.size() != Y.size() ){
    Rcpp::stop("Problem with unequal vector sizes when doing linear interpolation.");}
  //cout <<  " X(0): " <<  X(0) <<" X(Y.size()-1): " <<X(Y.size()-1)   <<   " Point of interest: " << X_PointOfInterest<< endl;
  if ( X_PointOfInterest < X(0) || X_PointOfInterest > X(Y.size()-1) ){Rcpp::warning("You interpolate out of the curve boundaries"); return(-1);}
  //Find the points right before and right after the point of interest
  for (int i=1; i<X.size() ; i++){
    if (X(i)>= X_PointOfInterest){
      xkp1 = X(i);
      xk = X(i-1);
      ykp1 = Y(i);
      yk = Y(i-1);
      break;}
  }
  //point-slope form for a line formula
  float t = (X_PointOfInterest -xk)/(xkp1 -xk);
  float yPOI = (1-t) * yk + t * ykp1;  // estimate point of interest
  // cout << "(" << xk << ",  " << X_PointOfInterest << " , " << xkp1 << ") & (" << yk << ", " << yPOI  << ", " << ykp1 << ")"<< endl;
  return (yPOI);   
} 

// [[Rcpp::export]]
Eigen::VectorXd RcppPseudoApprox(  const Eigen::Map<Eigen::VectorXd> & X,  const Eigen::Map<Eigen::VectorXd> & Y,  const Eigen::Map<Eigen::VectorXd> & X_target){
  //evaluate Y_target for X_target given X and Y
  
  int N =  X_target.size();
  Eigen::VectorXd rr(N); 
  for (int i=0; i<N ; i++){  
    rr(i) =  LinearInterpolation(X, Y, X_target(i)) ;  }
  return(rr);
}
