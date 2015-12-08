#' Functional Principal Component Scores Plot
#'
#' This function will create using the first two FPC scores a set of convex hulls of the scores as these hulls are defined by the references.
#'
#' @param fpcaObj An FPCA class object returned by FPCA().
#' @param factor Inflation factor for the bag-plot defining the loop of bag-plot or multiplying factor the KDE pilot bandwidth matrix.  (see ?aplpack::compute.bagplot/?ks::Hpi; default: 2.58/2).
#' @param variant Strng defining the outlier method used ('KDE' or 'bagplot') (default: 'KDE')
#' @param outlierList logical speciifying if a list with the grouping of points should be return (default: FALSE)
#' @param ... Additional arguments for the 'plot' function.
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createOutliersPlot(res)
#' @references
#' \cite{P. J. Rousseeuw, I. Ruts, J. W. Tukey (1999): The bagplot: a bivariate boxplot, The American Statistician, vol. 53, no. 4, 382-387}
#' \cite{R. J. Hyndman and H. L. Shang. (2010) Rainbow plots, bagplots, and boxplots for functional data, Journal of Computational and Graphical Statistics, 19(1), 29-45}
#'
#' @export

createOutliersPlot <- function(fpcaObj,variant = 'KDE',  factor = NULL, outlierList = FALSE, ...){
  
  if( !any( variant == c('KDE','bagplot')) ){
    stop("You request an outlier detection method.")
  }
  if ( variant == 'bagplot' && !is.element('aplpack', installed.packages()[,1]) ){
    stop("Cannot the use the bagplot method; the package 'aplpack' is unavailable.")
  }
  if ( variant == 'KDE' && !is.element('ks', installed.packages()[,1]) ){
    stop("Cannot the use the KDE method; the package 'ks' is unavailable.")
  } 
  
  
  
  xedge = 1.05 * max( abs(fpcaObj$xiEst[,1]))
  yedge = 1.05 * max( abs(fpcaObj$xiEst[,2]))  
  
  args1 <- list(  pch=10,  xlab=paste('FPC1 scores ', round(fpcaObj$cumFVE[1]), '%', sep=''   ),
                  ylab=paste('FPC2 scores ', round( diff( fpcaObj$cumFVE[1:2])), '%', sep='' ),  
                  xlim = c(-xedge, xedge), ylim =c(-yedge, yedge), lwd= 2)
  inargs <- list(...)
  args1[names(inargs)] <- inargs 
  
  if(length(fpcaObj$lambda) <2 ){
    stop("This plotting utility function needs at least two eigenfunctions.")
    return(NULL)
  } 
  
  if ( variant == 'bagplot' ){
    
    if ( is.null((factor))){
      factor = 2.58
    } 
    bgObj = aplpack::compute.bagplot(x= fpcaObj$xiEst[,1], y= fpcaObj$xiEst[,2], 
                                     approx.limit=3333 , factor = factor)     
    
    args2 = list (x = fpcaObj$xiEst[,1], y = fpcaObj$xiEst[,2], cex= .33,  panel.first = grid(), type='n' )
    
    do.call(plot, c(args2, args1))   
    # I do this because panel.first() does not work with the do.call()
    points(x = fpcaObj$xiEst[,1], y = fpcaObj$xiEst[,2], cex= .33,  panel.first = grid(),  lwd= 2) 
    
    lines( bgObj$hull.bag[c(1:nrow(bgObj$hull.bag),1),], col=2, lwd=2)
    lines( bgObj$hull.loop[c(1:nrow(bgObj$hull.loop),1),], col=4, lwd=2) 
    legend(legend= c('0.500', 'The fence'), x='topright', col=c(2,4), lwd=2)
    
    if(outlierList){
      return( list( 'bag' = match( apply(bgObj$pxy.bag,1, prod), apply( bgObj$xydata,1, prod)),
                    'loop'= match( apply(bgObj$pxy.outer,1, prod), apply( bgObj$xydata,1, prod)), 
                    'outlier' = match( apply(bgObj$pxy.outlier,1, prod) ,apply( bgObj$xydata,1, prod))) )
    } 
  } else {
    
    if ( is.null((factor))){
      factor = 2
    } 
    fhat <- ks::kde(x=fpcaObj$xiEst[,1:2], gridsize = c(400,400), compute.cont = TRUE, 
                    H = ks::Hpi( x=fpcaObj$xiEst[,1:2], binned=TRUE,  pilot="dscalar"  ) *  factor) 
    
    maxIndex = which( fhat$estimate == max(fhat$estimate), arr.ind = TRUE)
    qq = quickNNeval(xin = fhat$eval.points[[1]], yin = fhat$eval.points[[2]], 
                     # zin =  fhat$estimate,  
                     zin = monotoniseMatrix( fhat$estimate, maxIndex[1], maxIndex[2]), 
                     xout = fpcaObj$xiEst[,1], yout = fpcaObj$xiEst[,2] ) 
    
    args2 = list (x= fhat$eval.points[[1]], y=fhat$eval.points[[2]], z = monotoniseMatrix( fhat$estimate, NULL, NULL), 
                  labcex=1.66, col= c('black','blue','red'), levels = fhat$cont[c(50, 95, 99)], labels = c('50%', '95%', '99%'))
    do.call(contour, c(args2, args1)); 
    grid()
    
    points(fpcaObj$xiEst[qq <=  fhat$cont[99],1:2],cex=0.5, col='orange', pch=10 , lwd =2 ) 
    points(fpcaObj$xiEst[qq >  fhat$cont[99] & qq <=  fhat$cont[95],1:2],cex=0.33, col='red', pch=10, lwd =2 ) 
    points(fpcaObj$xiEst[qq >  fhat$cont[95] & qq <=  fhat$cont[50],1:2],cex=0.33, col='blue', pch=10 , lwd =2 ) 
    points(fpcaObj$xiEst[qq >=  fhat$cont[50],1:2],cex=0.33, col='black' , pch=10, lwd =2 )
    
    legend('bottomleft', c('< 50%','50%-95%','95%-99%','> 99%'), pch = 19, 
           col= c('black','blue','red', 'orange'), pt.cex=1.5, bg='white' )
    
    if (outlierList){
      return( list( 'p0to50'= which(qq >=  fhat$cont[50]),
                    'p50to95' = which(qq >  fhat$cont[95] & qq <=  fhat$cont[50]),
                    'p95to99' = which(qq >  fhat$cont[99] & qq <=  fhat$cont[95]),
                    'p99plus' = which(qq <=  fhat$cont[99]) ))
    }
  } 
  
}

quickNNeval <- function(xin,yin, zin, xout, yout){
  xindeces = sapply( xout, function(myArg) which.min( abs( xin - myArg) ), simplify = TRUE)
  yindeces = sapply( yout, function(myArg) which.min( abs( yin - myArg) ), simplify = TRUE )
  return( zin[ cbind(xindeces,yindeces)] )
}

monotonise <- function(x, maxIndex = NULL){
  xq = x;
  if (is.null(maxIndex)){
    maxIndex = which.max(x);
  }
  
  if( maxIndex != length(x) ){
    for (i in 1:( length(x) - maxIndex)){
      if( xq[ i + maxIndex] > xq[maxIndex + i - 1] ){
        xq[ i + maxIndex] = xq[maxIndex + i - 1]
      }
    }
  }
  if (maxIndex >= 3){
    for (i in 1:(maxIndex - 2 )){
      if( xq[ - 1 - i + maxIndex] > xq[maxIndex - i] ){
        xq[ - 1- i + maxIndex] = xq[maxIndex - i]
      }
    }
  }
  return(xq)
} 

monotoniseMatrix = function(zin, xmaxind, ymaxind){
  if(is.null(xmaxind) && is.null(ymaxind)){
    maxIndx = which( max(zin) == zin, arr.ind = TRUE)
    xmaxind = maxIndx[1]
    ymaxind = maxIndx[2]
  }
  zq = zin;
  for (j in 1:dim(zin)[2]){
    for (i in 1:dim(zin)[1]){
      if (i == 1 || j == 1 || j == dim(zin)[1] || i == dim(zin)[2]){
        sizeOut = max( abs(xmaxind - i) +1, abs(ymaxind - j) +1 )
        xcoord = round(   ( seq(i, xmaxind , length.out = sizeOut) ) )
        ycoord = round(   ( seq(j, ymaxind , length.out = sizeOut) ) ) 
        zq[ cbind(xcoord,ycoord) ] = monotonise( zq[ cbind(xcoord,ycoord) ]) 
      }
    }
  }
  return(zq)
}

