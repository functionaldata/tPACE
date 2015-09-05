createFuncBoxPlot <- function(fpcaObj, addIndx =NULL, openNewDev = TRUE, variant= 'bagplot'){
 
  if(openNewDev){ 
    dev.new(width=6.95, height=5.0, noRStudioGD=TRUE) ; 
  }

  fittedCurves <- fitted(fpcaObj)   
  s <- fpcaObj$workGrid
  N <- nrow(fittedCurves)
  
  if ( variant == 'bagplot' && !is.element('aplpack', installed.packages()[,1])){
    warning('Cannot use bagplot because aplpack::compute.bagplot is unavailable; reverting to point-wise');
    variant = 'pointwise'
  }
  
  if ( variant == 'bagplot' && is.element('aplpack', installed.packages()[,1]) ){
  
    plot(type='n', s, s, ylim=range(fittedCurves), xlab='s', ylab='y(s)')   
     
    grid()         
    bgObj = aplpack::compute.bagplot(x= fpcaObj$xiEst[,1], y= fpcaObj$xiEst[,2], approx.limit=3333)     
    fittedCurvesFence = fittedCurves[ is.element( rowSums(fpcaObj$xiEst[,1:2]), rowSums(bgObj$pxy.outer) ),]; 
    fittedCurvesBag = fittedCurves[ is.element( rowSums(fpcaObj$xiEst[,1:2]), rowSums(bgObj$pxy.bag) ),];
   
    polygon(x=c(s, rev(s)), y = c(apply(fittedCurvesFence,2, min), 
            rev(apply(fittedCurvesFence,2, max))), col= 'lightgrey',border=0)
    polygon(x=c(s, rev(s)), y = c(apply(fittedCurvesBag,2, min), 
            rev(apply(fittedCurvesBag,2,max))), col= 'darkgrey',border=1)  
    lines(x=s, y= apply(fittedCurves,2, mean) , col='red')
  } else if (variant== 'pointwise'){
    plot(type='n', s, s, ylim=range(fittedCurves), xlab='s', ylab='y(s)')  
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
