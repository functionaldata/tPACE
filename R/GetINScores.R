# This function obtains the FPC scores for dense
# regular functional data by trapezoidal rule integration (see 
# https://en.wikipedia.org/wiki/Functional_principal_component_analysis)

######
# Input:
######  
# yvec: length p vector of dense regular functional observations 
# tvec: length p vector of observed time grids for the functional observations
######
# Output: 
######
# ret:   xiEst: n by length(lambda) matrix of estimated FPC scores
#        fittedY: n by p matrix of fitted/recovered functional observations
##########################################################################

GetINScores <- function(yvec, tvec, optns,obsGrid, mu, lambda, phi, sigma2=NULL){
  
  if(is.vector(phi)){
    phi=matrix(as.numeric(phi),nrow=length(phi),ncol=1)
  }
  
  if(length(lambda) != ncol(phi)){
    stop('No. of eigenvalues is not the same as the no. of eigenfunctions.')
  }
  
  #tau = sort(unique(signif( unlist(t),14 ))) # get observed time grid
  ranget <- diff(range(tvec))
  mu= approx(obsGrid,mu,tvec)$y
  cy = yvec - mu
  phi = apply(phi,2,function(phivec){return(approx(obsGrid,phivec,tvec)$y)})
  
  if(!is.matrix(phi)){
    phi=matrix(as.numeric(phi),nrow=1,ncol=length(phi))
  }
  
  xiEst = matrix(0,length(lambda)) 
  # Get Scores xiEst
  for(i in 1:length(lambda)){
    temp = cy * phi[,i]
    xiEst[i,1] = trapzRcpp(X = tvec[!is.na(temp)], Y = temp[!is.na(temp)])
    if (optns[['shrink']] && !is.null(sigma2)) {
      xiEst[i,1] <- xiEst[i,1] * lambda[i] / 
        (lambda[i] + ranget * sigma2 / length(tvec))
    }
  }
  
  # Get Fitted Y: n by p matrix on observed time grid
  fittedY = mu + t(phi %*% xiEst)
  
  ret = list('xiEst' = xiEst,xiVar=matrix(NA, length(lambda), length(lambda)), 'fittedY' = fittedY)
  
  return(ret)
  
}
