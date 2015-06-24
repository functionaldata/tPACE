# This function obtains the FPC scores for dense
# regular functional data by trapezoidal rule integration

######
# Input:
######  
# ymat: n by p matrix of dense regular functional observations 
# t: list of observed time grids for the functional observations
######
# Output: 
######
# ret: a list of:
#        xiEst: n by length(lambda) matrix of estimated FPC scores
#        fittedY: n by p matrix of fitted/recovered functional observations
##########################################################################

GetINScores <- function(ymat, t, optns, mu, lambda, phi){
  if(length(lambda) != ncol(phi)){
    stop('No. of eigenvalues is not the same as the no. of eigenfunctions.')
  }

  n = nrow(ymat)
  t = sort(unique(unlist(t))) # get observed time grid
  mumat = matrix(rep(mu, n), nrow = n, byrow = TRUE)
  cymat = ymat - mumat

  xiEst = matrix(0, nrow = n, ncol = length(lambda))

  # Get Scores xiEst
  for(i in 1:length(lambda)){
  	tempmat = cymat * matrix(rep(phi[,i],n), nrow = n, byrow = TRUE)
    xiEst[,i] = sapply(1:n, function(j) trapz(x = t, y = tempmat[j,]))
  }

  # Get Fitted Y: n by p matrix on observed time grid
  fittedY = mumat + t(phi %*% t(xiEst))

  ret = list('xiEst' = xiEst, xiVar = NULL, 'fittedY' = fittedY)

  return(ret)
}