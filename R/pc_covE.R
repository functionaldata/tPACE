# This function outputs the smoothed covariance function 
# along the diagonal with and without measurement error,
# together with the estimated variance of measurement error

######
# Input: 
######
# obsGrid:    vector of all observed time/measurement points in increasing order
# regGrid:  vector of output time-point-grid
# bw_userCov:  2-d vector, bandwidths along 2 directions for covariance surface smoothing
# rotationCut:  2-element vector in [0,1] indicating the percent of data truncated during 
#               sigma^2 estimation (default c(1/4,3/4))
# kernel:  kernel function used for 2d smoothing, default is 'epan'
# rcov:    a struct/list from function GetRawCov

######
# Output: a list/struct of
######
# sigma2:  estimated variance of measurement error
# xvar:    smoothed cov along diagonal without measurement error
# yvar:   smoothed cov along diagonal with measurement error

pc_covE = function(obsGrid, regGrid, bw_userCov, rotationCut=c(0, 1), kernel = 'epan', rcov){
  a0 = min(obsGrid)
  b0 = max(obsGrid)
  lint = b0 - a0
  
  rcutprop = rotationCut[2] - rotationCut[1]
  if(rcutprop <= 0 || rcutprop > 1){
    warning("Invalid option: rotationCut.")
  }
  rcutGrid = regGrid[intersect(which(regGrid > a0 + lint * rotationCut[1]),
    which(regGrid < a0 + lint * rotationCut[2]))]
  #out22 = rcutGrid # Xiongtao: out22 is redundant? Pantelis: Probably yes.

  tPairs = rcov$tPairs # time points pairs for raw covariance
  rcovdiag = rcov$diag # get raw covariance along diagonal direction

  # get smoothed covariance surface for x(t) using lwls2d

  cxxn = rcov$cxxn # off-diagonal terms

  if(length(rcov$count) != 0 && class(rcov)[1] != 'BinnedRawCov'){
    # for dataType="RegularwithMV" case, the raw covariance
    # matrix needs to be divided by the number of 
    # individual sums for each element in the matrix.
    # for dataType="Dense" case, the divider is n for
    # each subject.
    cxxn = cxxn / rcov$count
  }

  if (class(rcov)[1] == 'BinnedRawCov')
    win1 <- rcov$count
  else    
    win1 = rep(1, length(cxxn))

  # get smoothed variance function for y(t) (observed) using lwls1d
  if (class(rcov)[1] == 'BinnedRawCov') {
    rcovdiag <- cbind(rcov$tDiag, rcov$diagMeans)
    win2 <- rcov$diagCount
    cxxn <- rcov$meanVals
  } else
  win2 = rep(1, nrow(rcovdiag))

  # yvar is the smoothed variance function along the diagonal line
  # yvar = lwls1d(bw = bw_userCov[1], kern = kernel, xin = rcovdiag[,1],
  #  yin = rcovdiag[,2], win = win2, xout = rcutGrid, returnFit = FALSE)
  xorder = order(rcovdiag[,1]); 
    yvar = lwls1d(bw = bw_userCov[1], kernel_type = kernel, xin = rcovdiag[xorder,1],
                  yin = rcovdiag[xorder,2], win = win2, xout = rcutGrid)




  # Estimate variance of measurement error term
  # use quadratic form on diagonal to estimate Var(x(t))
  # xvar = rotateLwls2dV2(bw = bw_userCov[1], kern = kernel, xin = tPairs, yin = cxxn, win = win1, xout = cbind(rcutGrid, out22))
  xvar = rotateLwls2dV2(bw = bw_userCov[1], kern = kernel, xin = tPairs, yin = cxxn, win = win1, xout =  rcutGrid)
  sigma2 = trapzRcpp(rcutGrid, (yvar - xvar)) / (lint * rcutprop)

  if(sigma2 < 0){
    warning("Warning: estimated sigma2 is negative, reset to 1e-6 now!")
    sigma2 = 1e-6
  }

  return(list('sigma2' = sigma2, 'xvar' = xvar, 'yvar' = yvar))
}
