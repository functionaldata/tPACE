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
#  rho: regularization parameter for sigma2
#  fitLambda: eigenvalues by least squares fit method
######
# Output: 
######
#  ret: FPCA output object
##########################################################################

MakeResultFPCA <- function(optns, smcObj, mu, scsObj, eigObj, 
                           scoresObj, obsGrid, workGrid, rho=NULL, fitLambda=NULL, inputData){
  if(optns$dataType == 'Sparse'){
    ret <- list(sigma2 = scsObj$sigma2, lambda = eigObj$lambda, phi = eigObj$phi,
                xiEst = t(do.call(cbind, scoresObj[1, ])), xiVar = scoresObj[2, ], 
                # fittedY = scoresObj[3, ], 
                obsGrid = obsGrid, mu = mu, workGrid = workGrid, smoothedCov = scsObj$smoothCov, FVE = eigObj$cumFVE[eigObj$kChoosen], 
                cumFVE =  eigObj$cumFVE, fittedCov = eigObj$fittedCov, optns = optns, bwMu = smcObj$bw_mu, bwCov = scsObj$bwCov, rho=rho, fitLambda=fitLambda)
  } else if(optns$dataType %in% c('Dense','DenseWithMV')){
    ret <- list(sigma2 = scsObj$sigma2, lambda = eigObj$lambda, phi = eigObj$phi,
                xiEst = scoresObj$xiEst, xiVar = scoresObj$xiVar, 
                # fittedY = scoresObj$fittedY, 
                obsGrid = obsGrid, mu = mu, workGrid = workGrid, smoothedCov = scsObj$smoothCov, FVE = eigObj$cumFVE[eigObj$kChoosen], 
                cumFVE = eigObj$cumFVE, fittedCov = eigObj$fittedCov, optns = optns, bwMu = smcObj$bw_mu, bwCov = scsObj$bwCov)
  } else {
    stop('Other dataType choices not implemented yet!')
  }
  
  #if (!is.null(optns$maxK)){
  #  if( optns$maxK < length(ret$lambda) ){
  #    Knew = optns$maxK;
  #    ret$lambda <- ret$lambda[1:Knew];
  #    ret$phi <-  ret$phi[,1:Knew];
  #    ret$xiEst <-  ret$xiEst[,1:Knew]; 
  #    ret$xiVar <- lapply( ret$xiVar, function(x) x[1:Knew, 1:Knew])
  #  } else {
  #    warning("The number of components requested is higher than estimated number of components. Consider increasing the FVE threshold.")
  #  }
  #}
  
  if(!optns$lean){
    ret$inputData <- inputData;
  } else {
    ret$inputData <- NULL
  }
  class(ret) <- 'FPCA'
  return(ret)
}
