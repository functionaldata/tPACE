# This function makes the output object for function FPCA

######
# Input:
######  
#  optns: options for FPCA function
#  smcObj: smooth mean curve object
#  scsObj: smooth cov surface object
#  eigObj: eigen analysis object
#  scoresObj: FPC scores object
#  obsGrid: observed time grid
#  workGrid: time grid for smoothed covariance surface
######
# Output: 
######
#  ret: FPCA output object
##########################################################################

MakeResultFPCA <- function(optns, smcObj, mu, scsObj, eigObj,
	scoresObj, obsGrid, workGrid, rho=NULL){
  if(optns$dataType == 'Sparse'){
  	ret <- list(sigma2 = scsObj$sigma2, lambda = eigObj$lambda, phi = eigObj$phi,
  	  xiEst = t(do.call(cbind, scoresObj[1, ])), xiVar = scoresObj[2, ], 
      # fittedY = scoresObj[3, ], 
      obsGrid = obsGrid, mu = mu, workGrid = workGrid, smoothedCov = scsObj$smoothCov, 
      fittedCov = eigObj$fittedCov, optns = optns, bwMu = smcObj$bw_mu, bwCov = scsObj$bwCov, rho=rho)
  } else if(optns$dataType == 'Dense'){
  	ret <- list(sigma2 = scsObj$sigma2, lambda = eigObj$lambda, phi = eigObj$phi,
  	  xiEst = scoresObj$xiEst, xiVar = scoresObj$xiVar, 
      # fittedY = scoresObj$fittedY, 
      obsGrid = obsGrid, mu = mu, workGrid = workGrid, smoothedCov = scsObj$smoothCov, 
      fittedCov = eigObj$fittedCov, optns = optns, bwMu = smcObj$bw_mu, bwCov = scsObj$bwCov)
  } else {
  	stop('Other dataType choices not implemented yet!')
  }
  return(ret)
}