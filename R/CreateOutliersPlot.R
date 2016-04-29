#' Functional Principal Component Scores Plot using 'bagplot' or 'KDE' methodology
#'
#' This function will create, using the first two FPC scores, a set of convex hulls of the scores based on 'bagplot' or 'KDE' methodology.
#'
#' @param fpcaObj An FPCA class object returned by FPCA().   
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @param ... Additional arguments for the 'plot' function. 
#'
#' @details Available control options are 
#' \describe{
#' \item{ifactor}{inflation ifactor for the bag-plot defining the loop of bag-plot or multiplying ifactor the KDE pilot bandwidth matrix. (see ?aplpack::compute.bagplot; ?ks::Hpi respectively; default: 2.58; 2 respectively).}
#' \item{variant}{string defining the outlier method used ('KDE', 'NN' or 'bagplot') (default: 'KDE')}
#' \item{unimodal}{logical specifying if the KDE estimate should be unimodal (default: FALSE, relavant only for variant='KDE')}
#' \item{nSlices}{integer between 3 and 16, denoting the number of slices to be used (default: 4, relavant only for groupingType='slice') }
#' \item{colSpectrum}{character vector to be use as input in the 'colorRampPalette' function defining the outliers colous(default: c("red",  "yellow", 'blue'), relavant only for groupingType='slice') }
#' \item{groupingType}{string specifying if a slice grouping ('slice') or a standard percentile/bagplot grouping ('standard') should be returned (default: 'standard')} 
#' }
#'
#' @return An (temporarily) invisible copy of a list containing the labels associated with each of sample curves. 
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' CreateOutliersPlot(res)
#' @references
#' \cite{P. J. Rousseeuw, I. Ruts, J. W. Tukey (1999): The bagplot: a bivariate boxplot, The American Statistician, vol. 53, no. 4, 382-387}
#' \cite{R. J. Hyndman and H. L. Shang. (2010) Rainbow plots, bagplots, and boxplots for functional data, Journal of Computational and Graphical Statistics, 19(1), 29-45}
#'
#' @export

CreateOutliersPlot <- function(fpcaObj, optns = NULL, ...){
  
  
  if(is.null(optns$ifactor)){
    ifactor = NULL
  } else {
    ifactor = optns$ifactor
  }
  
  if(is.null(optns$outlierList)){
    outlierList = NULL
  } else {
    outlierList = optns$outlierList
  }
  
  if(is.null(optns$unimodal)){
    unimodal = NULL
  } else {
    unimodal = optns$unimodal
  }
  
  if(is.null(optns$colSpectrum)){
    colSpectrum = NULL
  } else {
    colSpectrum = optns$colSpectrum
  }
  
  if(is.null(optns$groupingType)){
    groupingType = 'standard'
  } else {
    groupingType = optns$groupingType
  }
  
  if(is.null(optns$variant)){
    variant = 'KDE'
  } else {
    variant = optns$variant
  }
  
  if(is.null(optns$nSlices)){
    nSlices = 4
  } else {
    nSlices = optns$nSlices
  }
  
  if( !any( groupingType == c('standard','slice')) ){
    stop("You request an groupingType method not currenty available.")
  }
  if( !any( variant == c('KDE','bagplot', 'NN')) ){
    stop("You request an outlier detection method not currenty available.")
  }
  if ( variant == 'bagplot' && !is.element('aplpack', installed.packages()[,1]) ){
    stop("Cannot the use the bagplot method; the package 'aplpack' is unavailable.")
  }
  if ( variant == 'KDE' && !is.element('ks', installed.packages()[,1]) ){
    stop("Cannot the use the KDE method; the package 'ks' is unavailable.")
  } 
  if ( !is.null(unimodal) && !is.logical(unimodal) ){
    stop("The variable 'unimodal' must be logical.")
  } 
  if (is.null(colSpectrum)){
    colFunc = colorRampPalette(c("red",  "yellow", 'blue'));
  } else {
    colFunc = colorRampPalette(colSpectrum)
  }
  if (!is.null(ifactor) && (1 >= ifactor) ){
    warning("It is nonsensical for an inflation factor to be <= 1. 'ifactor' set to 1.1.")
    ifactor = 1.1;
  }
  if ( !(3 <= nSlices) || !(16 >= nSlices) || !(nSlices %% 1 == 0) ){
    warning("nSlices must be between a natural number betweeb 3 and 16. 'nSlices' set to 4.")
    nSlices = 4;
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
    
    if ( is.null((ifactor))){ ifactor = 2.58 } 
    bgObj = aplpack::compute.bagplot(x= fpcaObj$xiEst[,1], y= fpcaObj$xiEst[,2], 
                                     approx.limit=3333 , factor = ifactor)     
    # I do this because panel.first() does not work with the do.call()
    
    if(groupingType =='standard'){
      
      args2 = list (x = fpcaObj$xiEst[,1], y = fpcaObj$xiEst[,2], cex= .33,  type='n' )
      do.call(plot, c(args2, args1))   
      points(x = fpcaObj$xiEst[,1], y = fpcaObj$xiEst[,2], cex= .33,  panel.first = grid(),  lwd= 2) 
      lines( bgObj$hull.bag[c(1:nrow(bgObj$hull.bag),1),], col=2, lwd=2)
      lines( bgObj$hull.loop[c(1:nrow(bgObj$hull.loop),1),], col=4, lwd=2) 
      legend(legend= c('0.500', 'The fence'), x='topright', col=c(2,4), lwd=2)
      
      return( invisible( list( 
        'bag' = match( apply(bgObj$pxy.bag,1, prod), apply( bgObj$xydata,1, prod)),
        'loop'= match( apply(bgObj$pxy.outer,1, prod), apply( bgObj$xydata,1, prod)), 
        'outlier' = ifelse( is.null(bgObj$pxy.outlier), NA, 
                            match( apply(bgObj$pxy.outlier,1, prod) ,apply( bgObj$xydata,1, prod)))
      ) ) ) 
    } else { # groupingType : slice
      
      N <- nrow(fpcaObj$xiEst[,1:2]) 
      kNNindeces95plus <- (1:N %in% match( apply(bgObj$pxy.outlier,1, prod) ,apply( bgObj$xydata,1, prod)))
      return( makeSlicePlot(nSlices, colFunc, p95plusInd = kNNindeces95plus, N, args1, 
                            xiEsts = fpcaObj$xiEst[,1:2] , lambdas = fpcaObj$lambda[1:2] ) )
      
    } 
  } else if (variant == 'KDE') {  # variant 'kde'
    
    if ( is.null((ifactor))){ ifactor = 2 } 
    fhat <- ks::kde(x=fpcaObj$xiEst[ ,1:2], gridsize = c(400,400), compute.cont = TRUE, 
                    H = ks::Hpi( x=fpcaObj$xiEst[ ,1:2], binned=TRUE,  pilot="dscalar"  ) *  ifactor) 
    zin = fhat$estimate
    if( !is.null(unimodal) && unimodal ){
      maxIndex = which( fhat$estimate == max(fhat$estimate), arr.ind = TRUE)
      zin = monotoniseMatrix( fhat$estimate, maxIndex[1], maxIndex[2])
    }
    qq = quickNNeval(xin = fhat$eval.points[[1]], yin = fhat$eval.points[[2]], zin =  zin, 
                     xout = fpcaObj$xiEst[,1], yout = fpcaObj$xiEst[,2] ) 
    
    if(groupingType =='standard'){ 
      args2 = list (x= fhat$eval.points[[1]], y=fhat$eval.points[[2]], z = zin, 
                    labcex=1.66, col= c('black','blue','red'), levels = fhat$cont[c(50, 95, 99)], labels = c('50%', '95%', '99%')) 
      do.call(contour, c(args2, args1)); 
      grid()
      
      points(fpcaObj$xiEst[qq <=  fhat$cont[99],1:2],cex=0.5, col='orange', pch=10 , lwd =2 ) 
      points(fpcaObj$xiEst[qq >  fhat$cont[99] & qq <=  fhat$cont[95],1:2],cex=0.33, col='red', pch=10, lwd =2 ) 
      points(fpcaObj$xiEst[qq >  fhat$cont[95] & qq <=  fhat$cont[50],1:2],cex=0.33, col='blue', pch=10 , lwd =2 ) 
      points(fpcaObj$xiEst[qq >=  fhat$cont[50],1:2],cex=0.33, col='black' , pch=10, lwd =2 ) 
      legend('bottomleft', c('< 50%','50%-95%','95%-99%','> 99%'), pch = 19, 
             col= c('black','blue','red', 'orange'), pt.cex=1.5, bg='white' )
      
      return( invisible( list( 'p0to50'= which(qq >=  fhat$cont[50]),
                               'p50to95' = which(qq >  fhat$cont[95] & qq <=  fhat$cont[50]),
                               'p95to99' = which(qq >  fhat$cont[99] & qq <=  fhat$cont[95]),
                               'p99plus' = which(qq <=  fhat$cont[99]) )))
    } else { # groupingType : slice 
      
      kNNindeces95plus <- qq <=  fhat$cont[95]
      return( makeSlicePlot(nSlices, colFunc, p95plusInd = kNNindeces95plus, N, args1, 
                            xiEsts = fpcaObj$xiEst[,1:2] , lambdas = fpcaObj$lambda[1:2] ) )
      
    }
  } else if (variant == 'NN') {
    
    centrePoint = c(0,0);
    distName = 'euclidean';
    N <- nrow(fpcaObj$xiEst[,1:2])
    k99 <- floor(0.99*N);
    k95 <- floor(0.95*N);
    k50 <- floor(0.50*N);
    scaledXi <- apply(fpcaObj$xiEst[,1:2], 2, scale)
    distances <- apply(scaledXi, 1, function(aRow) dist(x = rbind(aRow, centrePoint), method = distName) ) 
    kNNindeces0to99  <- sort(x = distances, index.return = TRUE)$ix[1:k99] # Partial sort should be better
    kNNindeces0to50  <- kNNindeces0to99[1:k50] 
    kNNindeces50to95 <- kNNindeces0to99[(1+k50):k95]  
    kNNindeces95to99 <- kNNindeces0to99[(1+k95):k99]  
    kNNindeces99plus <- setdiff(1:N, kNNindeces0to99)
    
    if(groupingType =='standard'){ 
      
      args2 = list (x = fpcaObj$xiEst[,1], y = fpcaObj$xiEst[,2], cex= .33,  type='n' )
      do.call(plot, c(args2, args1))   
      grid()
      points(fpcaObj$xiEst[kNNindeces99plus,1:2],cex=0.5, col='orange', pch=10 , lwd =2 ) 
      points(fpcaObj$xiEst[kNNindeces95to99,1:2],cex=0.33, col='red', pch=10, lwd =2 ) 
      points(fpcaObj$xiEst[kNNindeces50to95,1:2],cex=0.33, col='blue', pch=10 , lwd =2 ) 
      points(fpcaObj$xiEst[kNNindeces0to50, 1:2],cex=0.33, col='black' , pch=10, lwd =2 ) 
      legend('bottomleft', c('< 50%','50%-95%','95%-99%','> 99%'), pch = 19, 
             col= c('black','blue','red', 'orange'), pt.cex=1.5, bg='white' )
      
      return( invisible( list( 'p0to50'= kNNindeces0to50,
                               'p50to95' = kNNindeces50to95,
                               'p95to99' = kNNindeces95to99,
                               'p99plus' = kNNindeces99plus)))
    } else { # groupingType : slice
      
      kNNindeces95plus <- (1:N %in% setdiff(1:N, kNNindeces0to99[1:k95]))
      return( makeSlicePlot(nSlices, colFunc, p95plusInd = kNNindeces95plus, N, args1, 
                            xiEsts = fpcaObj$xiEst[,1:2] , lambdas = fpcaObj$lambda[1:2] ) )
      
    }
  }
}

makeSlicePlot <- function( nSlices, colFunc, p95plusInd, N, args1, args2, xiEsts, lambdas ){
  
  kNNindeces95plus <- p95plusInd
  
  args2 = list (x = xiEsts[,1], y = xiEsts[,2], cex= .33,  type='n' )
  do.call(plot, c(args2, args1))   
  grid()
  
  points(xiEsts[!p95plusInd,1:2],cex=0.33, col='black' , pch=10, lwd =2 )
  Qstr =xiEsts[,1:2]  /  matrix( c( rep( sqrt(lambdas[1:2]), 
                                         each= length(xiEsts[,1]))), ncol=2); 
  
  colPal = colFunc( nSlices )
  v = 1:nSlices;
  colPal = colPal[v] # this just gives a smooth change and maximized the diffference between oppositve slices
  outlierList <- list()
  steps =  seq(-1, (nSlices-1) *2 -1 , by =2 )
  for( i in 1:nSlices){
    angle =steps[i] * pi/nSlices
    multiplier1 = sign( sin( angle + pi/2) )
    multiplier2 = sign( cos( angle + pi/ (nSlices/2)))
    qrtIndx =  multiplier1 * Qstr[,2] > multiplier1 * tan(angle) * Qstr[,1] & 
      multiplier2 * Qstr[,2] < multiplier2 * tan(angle + pi/ (nSlices/2) ) * Qstr[,1] 
    outlierList[[i]] = qrtIndx & kNNindeces95plus
    points(xiEsts[ outlierList[[i]],1:2], cex=0.93, col= colPal[i], pch=3, lwd =2 )  
  }
  return( invisible( list(  'p0to95'= which(p95plusInd), 
                            'outlier' = sapply(outlierList, which),
                            'outlierColours' = colPal)) ) 
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
