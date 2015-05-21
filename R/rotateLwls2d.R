# Rotate the data and then smooth the diagonal elements. We use quadratic terms on either direction, rather than only orthogonal to the diagonal.
# xout: a matrix of two columns containing the diagonal elements.

rotateLwls2d <- function(bw, kern='epan', xin, yin, win=NULL, xout) {
    xin <- rotate(as.matrix(xin), pi / 4)
    xout <- rotate(as.matrix(xout), pi / 4)
    
    datin <- data.frame(cbind(xin, as.numeric(yin)))
    names(datin) <- c('x1', 'x2', 'y')
    
    if (missing(win) || is.null(win))
        win <- rep(1, nrow(datin))

    # decide space allocation
    # uniqPts <- sort(unique(datin$x1))
    # r <- diff(range(uniqPts))
    # maxk <- round((2 * bw / r * length(uniqPts))^2 * 4)
    
    # browser()
    # fit <- locfit(y ~ lp(x1, x2, h=bw, deg=2, scale=FALSE), data=datin, weights=win, kern=kern) # , maxk=maxk)
    # val <- predict(fit, xout)
    val <- smooth.lf(lp(xin[, 1], xin[, 2], h=bw), yin, xout, direct=FALSE, deg=2, maxk=500, kern=kern, weights=win)$y
    
    return(val)
}


# rotate points counter clockwise.
# mat: a matrix whose rows contain the points.
# rad: radius to rotate.
rotate <- function(mat, rad) {
    rmat <- matrix(c(cos(rad), sin(rad), -sin(rad), cos(rad)), 2, 2)
    return(mat %*% t(rmat))
}
