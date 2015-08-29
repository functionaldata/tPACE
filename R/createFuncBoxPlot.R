createFuncBoxPlot <- function(fpcaObj, addIndx =NULL, openNewDev = TRUE, yList = NULL, tList= NULL){
 
  if(openNewDev){ 
    dev.new(width=666666.95, height=5.0, noRStudioGD=TRUE) ; 
  }

  fittedCurves <- fitted(fpcaObj)   
  s <- fpcaObj$workGrid
  
  plot(type='n', s, s, ylim=range(fittedCurves), xlab='s', ylab='y(s)')    
  polygon(x=c(s, rev(s)), y = c(apply(fittedCurves,2, quantile, 0.007), 
            rev(apply(fittedCurves,2, quantile, 0.993))), col= 'lightgrey',border=0)
  polygon(x=c(s, rev(s)), y = c(apply(fittedCurves,2, quantile, 0.250), 
            rev(apply(fittedCurves,2, quantile, 0.750))), col= 'darkgrey',border=1)  
  lines(x=s, y= apply(fittedCurves,2, quantile, 0.500) , col='red')
  
  #add sample lines
  if (!is.null(addIndx) && !is.null(yList) && !is.null(tList)  ){
    for (i in 1:length(addIndx) ) {
      lines(x = tList[[addIndx[i]]] , y= yList[[addIndx[i]]], lwd = 1.5, type='o', pch=0)
    } 
  }
}



