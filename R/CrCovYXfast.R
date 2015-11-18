#' Functional Cross Covariance between longitudinal variable Y and longitudinal variable X
#' 
#' Calculate the raw and the smoothed cross-covariance between functional predictors using bandwidth bw or estimate that bw using GCV. 
#' 
#' @param fpcaX Object of FPCA class 
#' @param fpcaY Object of FPCA class 
#' @examples
#' Ly1= list( rep(2.1,7), rep(2.1,3),2.1 );
#' Lt1 = list(1:7,1:3, 1);
#' Ly2 = list( rep(1.1,7), rep(1.1,3),1.1); 
#' Lt2 = list(1:7,1:3, 1);
#' Ymu1 = rep(55,7);
#' Ymu2 = rep(1.1,7);
#' AA<-CrCovYX(Ly1 = Ly1, Ly2= Ly2, Lt1=Lt1, Lt2=Lt2, Ymu1=Ymu1, Ymu2=Ymu2)
#'   
#' @references
#' \cite{Yang, Wenjing, Hans‐Georg Müller, and Ulrich Stadtmüller. "Functional singular component analysis." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 73.3 (2011): 303-324}
#' @export

CrCovYXfast <- function( fpcaX, fpcaY ){
  
  KsiCov = cov( fpcaX$xiEst, fpcaY$xiEst); 
  CrCov  = matrix( rep(0, length(fpcaX$phi[,1]) * length(fpcaY$phi[,1]) ), length(fpcaX$phi[,1]))
  
  CorrectionTermX = fpcaX$lambda[1] / apply(fpcaX$xiEst,2,var)[1] 
  CorrectionTermY = fpcaY$lambda[1] / apply(fpcaY$xiEst,2,var)[1] 
   
  for(i in 1:nrow(KsiCov)){  
    for(j in 1:ncol(KsiCov)){  
      CrCov = CrCov +  
        CorrectionTermX * CorrectionTermY * KsiCov[i,j] * (fpcaX$phi[,i] %*%  t(fpcaY$phi[,j])) 
    }
  }
  
  abc1 = 4;
  
  return(CrCov) 
}
