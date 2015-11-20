#' Functional Principal Component Analysis mode of variation plot
#' 
#' This function will open a new device if not instructed otherwise.
#'
#' @param ret An FPCA class object returned by FPCA(). 
#' @param k The k-th mode of variation to plot (default k = 1)
#' @param titleString a string variable to be used as title
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createDiagnosticsPlot(res)
#' @export

createModeOfVarPlot <-function(ret,  k = 1, titleString = NULL){ 
  
  if(k> length(ret$lambda) ){
    stop("You are asking to plot a mode of variation that is incomputable.")
  }  
  
  obsGrid = ret$obsGrid      
  s = ret$workGrid
  mu = approx( y = ret$mu, x = obsGrid, xout = s)$y
  
  sigma = sqrt(ret$lambda[k])
  sigma1 = sqrt(ret$lambda[1])
  phi = ret$phi[,k]
  phi1 = ret$phi[,1]
  
  plot(type='n', s, s, ylim=range(c( 3* sigma1 * phi1 + mu , -3* sigma1 * phi1 + mu  )), xlab='s', ylab='y(s)', main = titleString)  
  grid()    
    
  polygon(x=c(s, rev(s)), y = c( -2* sigma * phi + mu, 
                                rev(2* sigma * phi + mu)), col= 'lightgrey',border=NA)
  polygon(x=c(s, rev(s)), y = c( -1* sigma * phi + mu, 
                                    rev(1* sigma * phi + mu)), col= 'darkgrey',border=NA)  
  lines(x=s, y=mu , col='red')
  
}

