# 1 dimensional local weighted least squares smoother.
# nwe argument (in the Matlab version) is merged into kern.
# bw: bandwidth, a scalar.
# kern: kernel used: 'gauss', 'rect', 'tria', 'epan' (default), 'quar', 'trwt', or 'tcub'. (See \url{http://en.wikipedia.org/wiki/Kernel_%28statistics%29#In_non-parametric_statistics})
# npoly: degree of local polynomial; 1 corresponds to local linear.
# nder: order of derivative.
# xin: a vector of x-coordinate.
# yin: a vector of y-coordinate.
# win: a vector of weights on the observations.
# xout: a vector of output x-coordinate grid.
# returnFit: If TRUE, return the fitted locfit object for the convenience of future operations such as gcv. Use predict(fit) to obtain the fitted values.
# Difference from lwls.m: There is no nwe argument. The kernels supported are different. No local bandwidth choice is supported as it is never actually called in Matlab PACE.
lwls1d <- function(bw, kern='epan', npoly=1, nder=0, xin, yin, win=NULL, xout=dat(), returnFit=FALSE) {
    if (npoly < nder)
        stop('Degree of Polynomial should be no less than the order of derivative')

    if (bw <= 0L)
        stop('Bandwidth choice must be positive')

    datin <- data.frame(x=as.numeric(xin), y=as.numeric(yin))
    
    if (nder == 0L)
        nder <- numeric(0) # as in locfit
    else if (nder == 2L)
        nder <- c(1, 1) # as in locfit

    if (missing(win))
        win <- rep(1, length(xin))

    fit <- locfit(y ~ lp(x, h=bw, deg=npoly), datin, weights=win, deriv=nder, ev=xout, kern=kern)

    if (returnFit)
        return(fit)

    if (missing(xout)) {
        ret <- fitted(fit)
    } else {
        ret <- predict(fit)
    }

    return(ret)
}

