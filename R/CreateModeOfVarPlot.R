#' Functional Principal Component Analysis mode of variation plot
#' 
#' Create the k-th mode of variation plot around the mean. The red-line is
#' the functional mean, the grey shaded areas show the range of variations
#' around the mean: \eqn{ \pm Q \sqrt{\lambda_k} \phi_k}{+/- Q sqrt{lambda_k} phi_k}
#' for the dark grey area Q = 1, and for the light grey are Q = 2. In the case of 'rainbowPlot'
#' the blue edge corresponds to Q = -3, the green edge to Q = +3 and the red-line to Q = 0 (the mean).
#'
#' @param fpcaObj An FPCA class object returned by FPCA(). 
#' @param k The k-th mode of variation to plot (default k = 1) 
#' @param rainbowPlot Indicator to create a rainbow-plot instead of a shaded plot (default: FALSE)
#' @param colSpectrum Character vector to be use as input in the 'colorRampPalette' function defining the outliers colours (default: c("blue","red", "green"), relavant only for rainbowPlot=TRUE)
#' @param ... Additional arguments for the \code{plot} function.
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' CreateModeOfVarPlot(res)
#' @export

CreateModeOfVarPlot <-function(fpcaObj,  k = 1, rainbowPlot = FALSE, colSpectrum = NULL, ...){ 
  
  args1 <- list( main="Default Title", xlab='s', ylab='')  
  inargs <- list(...)
  
  if(k> length(fpcaObj$lambda) ){
    stop("You are asking to plot a mode of variation that is incomputable.")
  }  
  
  if( is.null(colSpectrum) ){
    colSpectrum = c("blue","red", "green")
  } 
  
  
  
  obsGrid = fpcaObj$obsGrid      
  s = fpcaObj$workGrid
  mu = fpcaObj$mu
  
  sigma = sqrt(fpcaObj$lambda[k])
  sigma1 = sqrt(fpcaObj$lambda[1])
  phi = fpcaObj$phi[,k]
  phi1 = fpcaObj$phi[,1]
  
  if( is.null(inargs$ylim)){
    inargs$ylim=range(c( 3* sigma1 * phi1 + mu , -3* sigma1 * phi1 + mu ))
  }
  
  args1[names(inargs)] <- inargs
  
  if(!rainbowPlot){
    do.call(plot, c(list(type='n'), list(x=s), list(y=s), args1))
    grid()    
    polygon(x=c(s, rev(s)), y = c( -2* sigma * phi + mu, 
                                   rev(2* sigma * phi + mu)), col= 'lightgrey',border=NA)
    polygon(x=c(s, rev(s)), y = c( -1* sigma * phi + mu, 
                                   rev(1* sigma * phi + mu)), col= 'darkgrey',border=NA)  
    lines(x=s, y=mu , col='red') 
  } else {
    # The numver of modes as well as the colour palette could/shoud be possibly user-defined
    numOfModes = 257 # Just a large number of make things look "somewhat solid".
    thisColPalette = colorRampPalette(colSpectrum)(numOfModes)
    
    Qmatrix = (seq(-2 ,2 ,length.out = numOfModes) %*% t(fpcaObj$phi[,k] * sqrt(fpcaObj$lambda[k]))) + 
      matrix(rep(mu,numOfModes), nrow=numOfModes, byrow = TRUE)  
    do.call(matplot, c(list(type='l'), list(x=s), list(y=t(Qmatrix)), list(col= thisColPalette), list(lty=1), list(lwd=2), args1))
    grid()    
    lines(x=s, y=mu , col= median(thisColPalette)) 
  }
}
