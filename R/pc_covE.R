# This function outputs the smoothed covariance function 
# along the diagonal with and without measurement error,
# together with the estimated variance of measurement error

######
# Input: 
######
# out1:    vector of all observed time/measurement points in increasing order
# out21:  vector of output time-point-grid
# bw_xcov:  2-d vector, bandwidths along 2 directions for covariance surface smoothing
# cut:    do not cut (0) or cut (1) the domain on both boundaries for 
#      smoothing along the diagonal direction, default is 1
# kernel:  kernel function used for 2d smoothing, default is 'epan'
# rcov:    a struct/list from function GetRawCov

######
# Output: a list/struct of
######
# sigma2:  estimated variance of measurement error
# xvar:    smoothed cov along diagonal without measurement error
# yvar:   smoothed cov along diagonal with measurement error

pc_covE = function(out1, out21, bw_xcov, cut = 1, kernel = 'epan', rcov){
  a0 = min(out1)
  b0 = max(out1)
  lint = b0 - a0
  out22 = out21

  tpairn = rcov$tpairn # time points pairs for raw covariance
  rcovdiag = rcov$diag # get raw covariance along diagonal direction

  # get smoothed covariance surface for x(t) using lwls2d

  cxxn = rcov$cxxn # off-diagonal terms

  if(length(rcov$count) != 0){
    # for regular="RegularwithMV" case, the raw covariance
    # matrix needs to be divided by the number of 
    # individual sums for each element in the matrix.
    # for regular="Dense" case, the divider is n for
    # each subject.
    cxxn = cxxn / rcov.count
  }

  win1 = rep(1, length(cxxn))

  # get smoothed variance function for y(t) (observed) using lwls1d
  win2 = rep(1, nrow(rcovdiag))

  # yvar is the smoothed variance function along the diagonal line
  yvar = lwls1d(bw = bw_xcov[1], kern = kernel, xin = rcovdiag[,1],
    yin = rcovdiag[,2], win = win2, xout = out21, returnFit = FALSE)

  # Estimate variance of measurement error term
  # use quadratic form on diagonal to estimate Var(x(t))
  xvar = rotateLwls2d(bw = bw_xcov[1], kern = kernel, 
    xin = tpairn, yin = cxxn, win = win1, xout = cbind(out21, out22))

  #expgrid = expand.grid(xout1, xout2)
  #eqind1 = which(expgrid[,1] == expgrid[,2])

  if(cut == 0){
    sigma2 = trapz(out21, yvar - xvar) / lint
  } else if(cut == 1){
    a = a0 + lint * 0.25
    b = a0 + lint * 0.75
    ind1 = intersect(which(out21 > a), which(out21 < b))
    yvar1 = yvar[ind1]
    xvar1 = xvar[ind1]
    sigma2 = trapz(out21[ind1], yvar1 - xvar1) * 2 / lint
  }

  if(sigma2 < 0){
    warning("Warning: estimated sigma2 is negative, reset to zero now!")
    sigma2 = 0
  }

  return(list('sigma2' = sigma2, 'xvar' = xvar, 'yvar' = yvar))
}