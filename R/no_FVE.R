# Selects best number of principal components based on
# FVE (fraction of variance explained) for the time being
# Other criteria: AIC, or BIC will be implemented later

######
# Input:
######
# userCov:           nRegGrid * nRegGrid matrix of smooth covariance surface.
# FVEthreshold:  a positive number between 0 and 1, ndicating the number
#                 of selected PCs will explain at least this percentage of
#                 total variation. Default value: 0.85.
######
# Output:
######
# no_opt:         positive integer, the best number of PCs chosen with FVE criterion 
# FVE:            no_opt * 1 vector, the dth entry corresponds to the cumulative 
#                 percentage of variation explained by the first d PCs
# lambda:         no_opt * 1 vector of the positive eigenvalues in a decreasing order
# eVec:           (unnormalized) A nRegGrid * no_opt matrix of eigenvectors 
#                 obtained from the eigen decomposition of the smooth covariance
#                 matrix userCov
#
# Note that FVE and lambda might be smaller if some of the eigenvalues are negative
# or complex numbers.
# This function uses rARPACK library for eigen-decomposition

no_FVE <- function(userCov, FVEthreshold=0.9999, returnEVec=FALSE, verbose=FALSE){
  numGrids = nrow(userCov)
  #optns.v0 = t(seq(0.1,0.9,length.out = numGrids))
  eigObj <- eigen(userCov)
  d <- eigObj$values
  eVec <- eigObj$vectors
  
  if(any(is.complex(d))){ # if there are any imaginary eigenvalues
    stop(sprintf("Some eigenvalues are complex. The estimated auto-covairance surface is not symmetric!"))
  }
  
  idx <- (d <= 0) # to remove nonpositive eigenvalues
  if(sum(idx) > 0) { # if there are any nonpositive eigenvalues
    if (verbose)
      warning(sprintf("%d real eigenvalues are negative or zero and are removed!", sum(idx)))
  }
  
  d <- d[!idx] 
  eVec <- eVec[, !idx]
  #FVE <- cumsum(d) / sum(d) # cumulative FVE for all available eigenvalues from fitted cov
  #no_opt <- min(which(FVE >= FVEthreshold))

  lambda <- d
  
  if (!returnEVec)
    eVec <- NULL
    
  return(list(lambda=lambda, eVec=eVec))
}

