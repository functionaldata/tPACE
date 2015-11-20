#' Create functional boxplot using 'bagplot' or 'pointwise' methodology
#'
#' Using an FPCA object create a functional box-plot based on the function scores.
#'
#' @param fpcaObj A object of class FPCA returned by the function FPCA().
#' @param addInx A vector of indeces corresponding to which samples one should overlay (Default: NULL)
#' @param variant A character variable indicating which methodology should be used ('bagplot' or 'pointwise')to create the functional box-plot (Default: 'bagplot')
#' @param titleString a string variable to be used as title
#' 
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createFuncBoxPlot(res, addIndx=c(1:3)) 
#' @references
#' \cite{P. J. Rousseeuw, I. Ruts, J. W. Tukey (1999): The bagplot: a bivariate boxplot, The American Statistician, vol. 53, no. 4, 382-387}
#'
#' @export

createFuncBoxPlot <- function(fpcaObj, addIndx =NULL, variant= 'bagplot', titleString = NULL){
 
  #if(openNewDev){ 
  #  dev.new(width=6.95, height=5.0, noRStudioGD=TRUE) ; 
  #}

  if( is.na( any(match( variant, c('pointwise', 'bagplot') )) ) ){
   stop("This plotting utility function can only implement a 'bagplot' or 'pointwise' mapping.")
   return(NULL)
  }

  fittedCurves <- fitted(fpcaObj)   
  s <- fpcaObj$workGrid
  N <- nrow(fittedCurves)
  
  if ( variant == 'bagplot' && !is.element('aplpack', installed.packages()[,1])){
    warning('Cannot use bagplot because aplpack::compute.bagplot is unavailable; reverting to point-wise');
    variant = 'pointwise'
  }
 
  if ( length(fpcaObj$lambda) <2) {
    warning('There is a single component used.');
  }
 
  if ( variant == 'bagplot' && is.element('aplpack', installed.packages()[,1]) ){
  
    plot(type='n', s, s, ylim=range(fittedCurves, na.rm = TRUE), xlab='s', ylab='y(s)')     
    grid()         
    if (  length(fpcaObj$lambda) >1) {
      bgObj = aplpack::compute.bagplot(x= fpcaObj$xiEst[,1], y= fpcaObj$xiEst[,2], approx.limit=3333)     
      fittedCurvesFence = fittedCurves[ is.element( rowSums(fpcaObj$xiEst[,1:2]), rowSums(bgObj$pxy.outer) ),]; 
      fittedCurvesBag = fittedCurves[ is.element( rowSums(fpcaObj$xiEst[,1:2]), rowSums(bgObj$pxy.bag) ),];
    } else {
      bgObj = boxplot(plot=FALSE, fpcaObj$xiEst[,1] )
      fittedCurvesFence = fittedCurves[ (fpcaObj$xiEst > bgObj$stats[1]) & (fpcaObj$xiEst < bgObj$stats[5]),];
      fittedCurvesBag = fittedCurves[ (fpcaObj$xiEst > bgObj$stats[2]) & (fpcaObj$xiEst < bgObj$stats[4]),];
    }
    polygon(x=c(s, rev(s)), y = c(apply(fittedCurvesFence,2, min), 
            rev(apply(fittedCurvesFence,2, max))), col= 'lightgrey',border=0)
    polygon(x=c(s, rev(s)), y = c(apply(fittedCurvesBag,2, min), 
            rev(apply(fittedCurvesBag,2,max))), col= 'darkgrey',border=1)  
    lines(x=s, y= apply(fittedCurves,2, mean) , col='red')
  } else if (variant== 'pointwise'){
    plot(type='n', s, s, ylim=range(fittedCurves), xlab='s', ylab='y(s)', main = titleString)  
    grid()     
    polygon(x=c(s, rev(s)), y = c(apply(fittedCurves,2, quantile, 0.007), 
            rev(apply(fittedCurves,2, quantile, 0.993))), col= 'lightgrey',border=0)
    polygon(x=c(s, rev(s)), y = c(apply(fittedCurves,2, quantile, 0.250), 
            rev(apply(fittedCurves,2, quantile, 0.750))), col= 'darkgrey',border=1)  
    lines(x=s, y= apply(fittedCurves,2, quantile, 0.500) , col='red')
  } else  {
    stop('Additional variances are not yet implemented')
  }
 
  yList = fpcaObj$inputData$y
  tList = fpcaObj$inputData$t 
 
  #add sample lines
  if (!is.null(addIndx) && !is.null(yList) && !is.null(tList)  ){
    for (i in 1:length(addIndx) ) {
      lines(x = tList[[addIndx[i]]] , y= yList[[addIndx[i]]], lwd = 1.5, type='o', pch=0)
    } 
  }
}
