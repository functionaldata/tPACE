# Selects best number of principal components based on
# FVE (fraction of variance explained) for the time being
# Other criteria: AIC, or BIC will be implemented later

######
# Input:
######
# xcov:           ngrid * ngrid matrix of smooth covariance surface.
# FVE_threshold:  a positive number between 0 and 1, ndicating the number
#                 of selected PCs will explain at least this percentage of
#                 total variation. Default value: 0.85.
######
# Output:
######
# no_opt:         positive integer, the best number of PCs chosen with FVE criterion 
# FVE:            ngrid * 1 vector, the dth entry corresponds to the cumulative 
#                 percentage of variation explained by the first d PCs
# lambda:         ngrid * 1 vector of the positive eigenvalues in a decreasing order
#                 obtained from the eigen decomposition of the smooth covariance
#                 matrix xcov
#
# Note that FVE and lambda might be smaller if some of the eigenvalues are negative
# or complex numbers.
# This function uses rARPACK library for eigen-decomposition
library(rARPACK)
library(gtools) # for error and warning message

no_FVE <- function(xcov, FVE_threshold=0.85){
  numGrids = nrow(xcov)
  #opts.v0 = t(seq(0.1,0.9,length.out = numGrids))
  d = eigs(xcov, k = numGrids-2, which = "LR")$values
  # at most ngrid-2 eigenvalues can be obtained for nonsymmetric or complex problems
  # "LM" corresponds to Largest Magnitude, another option may be "LR": Largest Real part
  
  idx = which(is.complex(d)) # to remove any imaginary eigenvalues
  if(invalid(idx) == FALSE){ # if there are any imaginary eigenvalues
    stop(sprintf("%d eigenvalues are complex. The estimated auto-covairance surface is not symmetric!",length(idx)))
  }
  
  idx = which(d <= 0) # to remove nonpositive eigenvalues
  if(invalid(idx) == FALSE){ # if there are any nonpositive eigenvalues
    warning(sprintf("Warning: %d real eigenvalues are negative or zero and are removed!",length(idx)))
    d = d[d > 0]
  }
  
  lambda = d # eigenvalues to output
  FVE = cumsum(lambda)/sum(lambda) # cumulative FVE to output
  no_opt = min(which(FVE > FVE_threshold))
  return(list(no_opt, FVE, lambda))
}

