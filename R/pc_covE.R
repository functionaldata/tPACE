# This function outputs the smoothed covariance function 
# along the diagonal with and without measurement error,
# together with the estimated variance of measurement error

######
# Input: 
######
# t:		list of observed time points for n subjects (not used for now)
# out1:		vector of all observed time/measurement points in increasing order
# out21:	vector of output time-point-grid
# bw_xcov:	2-d vector, bandwidths along 2 directions for covariance surface smoothing
# cut:		do not cut (0) or cut (1) the domain on both boundaries for 
#			smoothing along the diagonal direction, default is 1
# kernel:	kernel function used for 2d smoothing, default is 'epan'
# rcov:		a struct/list from function GetRawCov

######
# Output: a list/struct of
######
# sigma:	estimated variance of measurement error
# xvar:		smoothed cov along diagonal without measurement error
# yvar: 	smoothed cov along diagonal with measurement error

library(caTools)
pc_covE = function(t, out1, out21, bw_xcov, cut = 1, kernel = 'epan', rcov){
	a0 = min(out1)
	b0 = max(out1)
	lint = b0 - a0
	out22 = out21

	tpairn = rcov$tpairn # time points pairs for raw covariance

	# get smoothed covariance surface for x(t) using lwls2d

	tneq = which(tpairn[1,] != tpairn[2,])
	cyy = rcov$cyy

	if(length(rcov$count) != 0){
		# for regular="RegularwithMV" case, the raw covariance
		# matrix needs to be divided by the number of 
		# individual sums for each element in the matrix.
		# for regular="Dense" case, the divider is n for
		# each subject.
		cyy = cyy / rcov.count
	}

	cxx = cyy[tneq] # off-diagnal terms
	win1 = rep(1, length(cxx))

	# get smoothed variance function for y(t) (observed) using lwls1d
	teq = which(tpairn[1,] == tpairn[2,])
	vyy = cyy[teq]
	win2 = rep(1, length(vyy))

	# yvar is the smoothed variance function along the diagonal line
	yvar = lwls1d(bw = bw_xcov[1], kern = kernel, xin = tpairn[1,teq],
		yin = vyy, win = win2, xout = out21, returnFit = FALSE)

	# Estimate variance of measurement error term
	# use quadratic form on diagonal to estimate Var(x(t))
	xvar = rotateLwls2d(bw = bw_xcov[1], kern = kernel, 
		xin = tpairn[,tneq], yin = cxx, win = win1, xout = cbind(out21, out22))

	#expgrid = expand.grid(xout1, xout2)
	#eqind1 = which(expgrid[,1] == expgrid[,2])

	if(cut == 0){
		sigma = trapz(out21, yvar - xvar) / lint
	} else if(cut == 1){
		a = a0 + lint * 0.25
		b = a0 + lint * 0.75
		ind1 = which(out21 > a && out21 < b)
		yvar1 = yvar[ind1]
		xvar1 = xvar[ind1]
		sigma = trapz(out21[ind1], yvar1 - xvar1) * 2 / lint
	}

	if(sigma < 0){
		warning("Warning: estimated sigma is negative, reset to zero now!")
		sigma = 0
	}

	return(list('sigma' = sigma, 'xvar' = xvar, 'yvar' = yvar))
}