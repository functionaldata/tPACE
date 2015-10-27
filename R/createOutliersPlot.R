#' Functional Principal Component Scores Plot
#'
#' This function will create using the first two FPC scores a set of convex hulls of the scores as these hulls are defined by the references.
#'
#' @param ret An FPCA class object returned by FPCA().
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createOutliersPlot(res)
#' @references
#' \item{P. J. Rousseeuw, I. Ruts, J. W. Tukey (1999): The bagplot: a bivariate boxplot, The American Statistician, vol. 53, no. 4, 382-387}
#'
#' @export

createOutliersPlot <- function(fpcaObj){
 
 # if( is.na( any(match( variant, c('pointwise', 'bagplot') )) ) ){
 #  stop("This plotting utility function can only implement a 'bagplot' or 'pointwise' mapping.")
 #  return(NULL)
 # }
  variant = 'bagplot'


  if(length(fpcaObj$lambda) <2 ){
    stop("This plotting utility function needs at least two eigenfunctions")
    return(NULL)
  }
  
   xedge = 1.1 * max( abs(fpcaObj$xiEst[,1]))
   yedge = 1.1 * max( abs(fpcaObj$xiEst[,2]))  
   xlabelString = paste('FPC1 scores ', round(fpcaObj$FVE[1]), '%', sep=''   )
   ylabelString = paste('FPC2 scores ', round( diff( fpcaObj$FVE[1:2])), '%', sep=''   )


  if ( 'obsolete' == 'pointwise'){
   
    plot(fpcaObj$xiEst[,1], fpcaObj$xiEst[,2], cex= .33, xlim = c(-xedge, xedge), ylim =c(-yedge, yedge), 
       pch=10, xlab= xlabelString, ylab=ylabelString )
    abline(v=(seq(-xedge, xedge, length.out=  21)), col="lightgray", lty="dotted")
    abline(h=(seq(-yedge, yedge, length.out=  21)), col="lightgray", lty="dotted")
    
    varXi1Xi2 = diag((fpcaObj$lambda[1:2]))
    muXi1Xi2 =  apply(cbind(fpcaObj$xiEst[,1], fpcaObj$xiEst[,2]),2,mean)

    lines(ellipse::ellipse( varXi1Xi2, centre= muXi1Xi2, level=.5), col=2)
    lines(ellipse::ellipse( varXi1Xi2, centre= muXi1Xi2, level=.95), col=3)
    lines(ellipse::ellipse( varXi1Xi2, centre= muXi1Xi2, level=.99), col=4)
    lines(ellipse::ellipse( varXi1Xi2, centre= muXi1Xi2, level=.999), col=5)   
    legend(legend= c('0.500', '0.950', '0.990', '0.999'), x='topright', col=2:5, lwd=2)
  } else if ( variant == 'bagplot' && is.element('aplpack', installed.packages()[,1]) ){
  
    bgObj = aplpack::compute.bagplot(x= fpcaObj$xiEst[,1], y= fpcaObj$xiEst[,2], approx.limit=3333)     
  
    plot(fpcaObj$xiEst[,1], fpcaObj$xiEst[,2], cex= .33, xlim = c(-xedge, xedge), ylim =c(-yedge, yedge), 
       pch=10, xlab= xlabelString, ylab=ylabelString )
    abline(v=(seq(-xedge, xedge, length.out=  21)), col="lightgray", lty="dotted")
    abline(h=(seq(-yedge, yedge, length.out=  21)), col="lightgray", lty="dotted")  
    lines( bgObj$hull.bag[c(1:nrow(bgObj$hull.bag),1),], col=2, lwd=2)
    lines( bgObj$hull.loop[c(1:nrow(bgObj$hull.loop),1),], col=4, lwd=2) 
    legend(legend= c('0.500', 'The fence'), x='topright', col=c(2,4), lwd=2)
  } else {
    stop('Additional outlier plots are not yet implemented')
  } 
 
}
