#' Functional Principal Component or Functional Singular Value Decomposition Scores Plot using 'bagplot' or 'KDE' methodology
#'
#' This function will create, using the first components scores, a set of convex hulls of the scores based on 'bagplot' or 'KDE' methodology.
#'
#' @param fObj A class object returned by FPCA() or FSVD().   
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @param ... Additional arguments for the 'plot' function. 
#'
#' @details Available control options are 
#' \describe{
#' \item{ifactor}{inflation ifactor for the bag-plot defining the loop of bag-plot or multiplying ifactor the KDE pilot bandwidth matrix. (see ?aplpack::compute.bagplot; ?ks::Hpi respectively; default: 2.58; 2 respectively).}
#' \item{variant}{string defining the outlier method used ('KDE', 'NN' or 'bagplot') (default: 'KDE')}
#' \item{unimodal}{logical specifying if the KDE estimate should be unimodal (default: FALSE, relevant only for variant='KDE')}
#' \item{maxVar}{logical specifying if during slicing we should used the directions of maximum variance (default: FALSE for FPCA, TRUE for FSVD)}
#' \item{nSlices}{integer between 3 and 16, denoting the number of slices to be used (default: 4, relevant only for groupingType='slice') }
#' \item{showSlices}{logical specifying if during slicing we should show the outline of the slice (default: FALSE)}
#' \item{colSpectrum}{character vector to be use as input in the 'colorRampPalette' function defining the outliers colours (default: c("red",  "yellow", 'blue'), relevant only for groupingType='slice') }
#' \item{groupingType}{string specifying if a slice grouping ('slice') or a standard percentile/bagplot grouping ('standard') should be returned (default: 'standard')} 
#' \item{fIndeces}{a two-component vector with the index of the mode of variation to consider (default: c(1,2) for FPCA and c(1,1) for FSVD)}
#' }
#'
#' @return An (temporarily) invisible copy of a list containing the labels associated with each of sample curves. 
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 420
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' CreateOutliersPlot(res)
#' }
#' @references
#' \cite{P. J. Rousseeuw, I. Ruts, J. W. Tukey (1999): The bagplot: a bivariate boxplot, The American Statistician, vol. 53, no. 4, 382-387}
#' \cite{R. J. Hyndman and H. L. Shang. (2010) Rainbow plots, bagplots, and boxplots for functional data, Journal of Computational and Graphical Statistics, 19(1), 29-45}
#' \cite{}
#' @export

CreateOutliersPlot <- function(fObj, optns = NULL, ...){
  
  fObjClass <- class(fObj)
  if( !(fObjClass %in% c('FSVD','FPCA')) ){
    stop("CreateOutliersPlot() expects an FPCA or an FSVD object as input.")
  } 
  
  newOptns <- CheckAndCreateCOPoptions(optns,fObjClass);
  
  nSlices = newOptns$nSlices;    ifactor = newOptns$ifactor; 
  colFunc = newOptns$colFunc;    fIndeces = newOptns$fIndeces; 
  variant = newOptns$variant;    groupingType = newOptns$groupingType; 
  unimodal = newOptns$unimodal;  outlierList = newOptns$outlierList;
  maxVar = newOptns$maxVar;      showSlices = newOptns$showSlices
  
  fVarAlls <- c();
  if(fObjClass == 'FPCA'){
    fVarAlls <-  fObj$lambda
  } else { 
    fVarAlls <-  (fObj$sValues)^2  
  }
  
  if(fIndeces[2] > length(fVarAlls)){
    stop("You requested a mode of variation that is not available.")
  }
  
  fScores1 <- fScores2 <- c();
  if(fObjClass == 'FPCA'){
    fScores1 <- fObj$xiEst[,fIndeces[1]]
    fScores2 <- fObj$xiEst[,fIndeces[2]]
  } else { 
    fScores1 <- fObj$sScores1[,fIndeces[1]]
    fScores2 <- fObj$sScores2[,fIndeces[2]]
  }
  fScoresAll <- cbind(fScores1, fScores2) 
  
  xedge = 1.05 * max( abs(fScores1))
  yedge = 1.05 * max( abs(fScores2))  
  
  args1 <- list();
  if(fObjClass == 'FSVD'){
    args1 <- list(  pch=10,  xlab=paste('S1 FSC', fIndeces[1] ,' scores ', sep=''   ),
                    ylab=paste('S2 FSC', fIndeces[2] ,' scores ', sep='' ),  
                    xlim = c(-xedge, xedge), ylim =c(-yedge, yedge), lwd= 2)
  } else {
    args1 <- list(  pch=10,  xlab=paste('FPC', fIndeces[1] ,' scores ', round(fObj$cumFVE[fIndeces[1]]), '%', sep=''   ),
                    ylab=paste('FPC', fIndeces[2] ,' scores ', round( diff( fObj$cumFVE[c(fIndeces[2]-1, fIndeces[2])])), '%', sep='' ),  
                    xlim = c(-xedge, xedge), ylim =c(-yedge, yedge), lwd= 2)
  } 
  inargs <- list(...)
  args1[names(inargs)] <- inargs
  
  nComp <- length(fVarAlls) 
  
  if(nComp <2 ){
    stop("This plotting utility function needs at least two functional components.")
    return(NULL)
  } 
  
  if ( variant == 'bagplot' ){
    
    if ( is.null((ifactor))){ ifactor = 2.58 } 
    bgObj = aplpack::compute.bagplot(x= fScores1, y= fScores2, approx.limit=3333 , factor = ifactor)     
    # I do this because panel.first() does not work with the do.call()
    
    if(groupingType =='standard'){
      
      args2 = list(x= fScores1, y= fScores2,  cex= .33,  type='n' )
      do.call(plot, c(args2, args1))   
      points(x = fScores1, y = fScores2, cex= .33,  panel.first = grid(),  lwd= 2) 
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
      
      N <- nrow(fScoresAll) 
      kNNindeces95plus <- (1:N %in% match( apply(bgObj$pxy.outlier,1, prod) ,apply( bgObj$xydata,1, prod)))
      return( makeSlicePlot(nSlices, colFunc, p95plusInd = kNNindeces95plus, N, args1, 
                            scoreEsts = fScoresAll , varEsts = fVarAlls[fIndeces], 
                            useDirOfMaxVar = maxVar, showSlices = showSlices) )
      
    } 
  } else if (variant == 'KDE') {  # variant 'kde'
    
    if ( is.null((ifactor))){ ifactor = 2 } 
    fhat <- ks::kde(x=fScoresAll, gridsize = c(400,400), compute.cont = TRUE, 
                    H = ks::Hpi( x=fScoresAll, binned=TRUE,  pilot="dscalar"  ) *  ifactor) 
    zin = fhat$estimate
    if( !is.null(unimodal) && unimodal ){
      maxIndex = which( fhat$estimate == max(fhat$estimate), arr.ind = TRUE)
      zin = monotoniseMatrix( fhat$estimate, maxIndex[1], maxIndex[2])
    }
    qq = quickNNeval(xin = fhat$eval.points[[1]], yin = fhat$eval.points[[2]], zin =  zin, 
                     xout = fScores1, yout = fScores2 ) 
    
    if(groupingType =='standard'){ 
      args2 = list (x= fhat$eval.points[[1]], y=fhat$eval.points[[2]], z = zin, 
                    labcex=1.66, col= c('black','blue','red'), levels = fhat$cont[c(50, 95, 99)], labels = c('50%', '95%', '99%')) 
      do.call(graphics::contour, c(args2, args1)); 
      grid(col = "#e6e6e6")
      
      points(fScoresAll[qq <=  fhat$cont[99], ],cex=0.5, col='orange', pch=10 , lwd =2 ) 
      points(fScoresAll[qq >  fhat$cont[99] & qq <=  fhat$cont[95], ],cex=0.33, col='red', pch=10, lwd =2 ) 
      points(fScoresAll[qq >  fhat$cont[95] & qq <=  fhat$cont[50], ],cex=0.33, col='blue', pch=10 , lwd =2 ) 
      points(fScoresAll[qq >=  fhat$cont[50], ],cex=0.33, col='black' , pch=10, lwd =2 ) 
      legend('bottomleft', c('< 50%','50%-95%','95%-99%','> 99%'), pch = 19, 
             col= c('black','blue','red', 'orange'), pt.cex=1.5, bg='white' )
      
      return( invisible( list( 'p0to50'= which(qq >=  fhat$cont[50]),
                               'p50to95' = which(qq >  fhat$cont[95] & qq <=  fhat$cont[50]),
                               'p95to99' = which(qq >  fhat$cont[99] & qq <=  fhat$cont[95]),
                               'p99plus' = which(qq <=  fhat$cont[99]) )))
    } else { # groupingType : slice 
      
      kNNindeces95plus <- qq <=  fhat$cont[95]
      return( makeSlicePlot(nSlices, colFunc, p95plusInd = kNNindeces95plus, N, args1, 
                            scoreEsts = fScoresAll , varEsts = fVarAlls[fIndeces], 
                            useDirOfMaxVar = maxVar, showSlices = showSlices) )
      
    }
  } else if (variant == 'NN') {
    
    centrePoint = c(0,0);
    distName = 'euclidean';
    N <- nrow(fScoresAll)
    k99 <- floor(0.99*N);
    k95 <- floor(0.95*N);
    k50 <- floor(0.50*N);
    scaledXi <- apply(fScoresAll, 2, scale)
    distances <- apply(scaledXi, 1, function(aRow) dist(x = rbind(aRow, centrePoint), method = distName) ) 
    kNNindeces0to99  <- sort(x = distances, index.return = TRUE)$ix[1:k99] # Partial sort should be better
    kNNindeces0to50  <- kNNindeces0to99[1:k50] 
    kNNindeces50to95 <- kNNindeces0to99[(1+k50):k95]  
    kNNindeces95to99 <- kNNindeces0to99[(1+k95):k99]  
    kNNindeces99plus <- setdiff(1:N, kNNindeces0to99)
    
    if(groupingType =='standard'){ 
      
      args2 = list (x = fScores1, y = fScores2, cex= .33,  type='n' )
      do.call(plot, c(args2, args1))   
      grid(col = "#e6e6e6")
      points(fScoresAll[kNNindeces99plus,],cex=0.5, col='orange', pch=10 , lwd =2 ) 
      points(fScoresAll[kNNindeces95to99,],cex=0.33, col='red', pch=10, lwd =2 ) 
      points(fScoresAll[kNNindeces50to95,],cex=0.33, col='blue', pch=10 , lwd =2 ) 
      points(fScoresAll[kNNindeces0to50, ],cex=0.33, col='black' , pch=10, lwd =2 ) 
      legend('bottomleft', c('< 50%','50%-95%','95%-99%','> 99%'), pch = 19, 
             col= c('black','blue','red', 'orange'), pt.cex=1.5, bg='white' )
      
      return( invisible( list( 'p0to50'= kNNindeces0to50,
                               'p50to95' = kNNindeces50to95,
                               'p95to99' = kNNindeces95to99,
                               'p99plus' = kNNindeces99plus)))
    } else { # groupingType : slice
      
      kNNindeces95plus <- (1:N %in% setdiff(1:N, kNNindeces0to99[1:k95]))
      return( makeSlicePlot(nSlices, colFunc, p95plusInd = kNNindeces95plus, N, args1, 
                            scoreEsts = fScoresAll, varEsts = fVarAlls[fIndeces],
                            useDirOfMaxVar = maxVar, showSlices = showSlices) )
      
    }
  }
}

makeSlicePlot <- function( nSlices, colFunc, p95plusInd, N, args1, args2, scoreEsts, varEsts, useDirOfMaxVar, showSlices){
  
  
  kNNindeces95plus <- p95plusInd
  
  args2 = list (x = scoreEsts[,1], y = scoreEsts[,2], cex= .33,  type='n' )
  do.call(plot, c(args2, args1))   
  grid(col = "#e6e6e6")
  
  points(scoreEsts[!p95plusInd, ],cex=0.33, col='black' , pch=10, lwd =2 )
  Qstr = apply(scoreEsts, 2, scale, center = FALSE) #
  #scoreEsts /  matrix( c( rep( sqrt(varEsts ), each= length(scoreEsts[,1]))), ncol=2); 
  
  dirOfMaxVar <- c(1,0);
  if(useDirOfMaxVar){
    dirOfMaxVar <- svd(scoreEsts, nv = 1)$v;
    if(all(dirOfMaxVar <0) ){
      dirOfMaxVar = -dirOfMaxVar
    }
    abline(0,  dirOfMaxVar[2]/dirOfMaxVar[1], col='magenta', lty=2) # Uncomment if you want to see the direction of max variance
  }
  
  colPal = colFunc( nSlices )
  v = 1:nSlices;
  colPal = colPal[v] # this just gives a smooth change and maximized the diffference between oppositve slices
  outlierList <- list()
  # steps =  seq(-1, (nSlices-1) *2 -1 , by =2 ) 
  angles <- seq(0,2*pi, length.out = nSlices + 1) -  1*pi/nSlices + atan2(dirOfMaxVar[2],dirOfMaxVar[1])
  sd1 = sd(scoreEsts[,1]); 
  sd2 = sd(scoreEsts[,2]); 
  for( i in 1:nSlices){
    angle = angles[i] # atan2(dirOfMaxVar[2],dirOfMaxVar[1]) + steps[i] * pi/nSlices
    multiplier1 = sign( sin( angle + pi/2) )
    multiplier2 = sign( cos( angle + pi/ (nSlices/2)))
    qrtIndx =  multiplier1 * Qstr[,2] > multiplier1 * tan(angle) * Qstr[,1] & 
      multiplier2 * Qstr[,2] < multiplier2 * tan(angle + pi/ (nSlices/2) ) * Qstr[,1] 
    outlierList[[i]] = qrtIndx & kNNindeces95plus
    points(scoreEsts[ outlierList[[i]], c(1,2), drop=FALSE], cex=0.93, col= colPal[i], pch=3, lwd =2 )
    if(showSlices){
      bigNumber = 10 * max(abs(as.vector(scoreEsts)))
      lines(x = c(0, bigNumber * multiplier1), col=colPal[i],  
            y = c(0, bigNumber * multiplier1 * tan(angle) * sd2 / sd1))
      #lines(x = c(0, bigNumber * multiplier2), col=colPal[i],  
      #      y = c(0, bigNumber * multiplier2 * tan(angle + pi/ (nSlices/2) ) * sd2 / sd1))
    }
  } 
  return( invisible( list(  'p0to95'= which(!p95plusInd), 
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
