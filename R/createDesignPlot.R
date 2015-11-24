#' Create the design plot of the functional data.
#'
#' This function will open a new device if not instructed otherwise.
#'
#' @param t a list of observed time points for functional data
#' @param obsGrid a vector of sorted observed time points
#' @param isColorPlot an option for colorful plot: 
#'                    TRUE: create color plot with color indicating counts
#'                    FALSE: create black and white plot with dots indicating observed time pairs
#' @param noDiagonal an option specifying plotting the diagonal design points:
#'                   TRUE:  remove diagonal time pairs
#'                   FALSE:  do not remove diagonal time pairs
#' @param yname the name of the variable containing functional observations
#' @param ... Other arguments passed into \code{plot()}. 
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' createDesignPlot(sampWiener$tList, sort(unique(unlist(sampWiener$tList))))
#' @export

createDesignPlot = function(t, obsGrid = NULL, isColorPlot=TRUE, noDiagonal=TRUE, yname = NULL, ...){
  if( is.null(obsGrid)){
    obsGrid = sort(unique(unlist(t)))
  }
  if( is.null(yname)){
    titleString =  paste('Design Plot')
  } else {
    titleString =  paste('Design Plot of', yname)
  }

  res = designPlotCount(t, obsGrid, noDiagonal, isColorPlot)
  # binrawcov is the function in place of designPlotCount
  # or getCount
  
  # construct plot structure
  #  plot(obsGrid[1], obsGrid[1], xlim = range(obsGrid) + isColorPlot * 0.25 * c(0,diff(range(obsGrid))),
  #   ylim = range(obsGrid),
  #   xlab = 'Observed time grid', ylab = 'Observed time grid',
  #   main = titleString, col = 'white')

  oldpty <- par()[['pty']]
  par(pty="s")
  if(isColorPlot == TRUE){
  	createColorPlot(res, obsGrid,titleString, ...)
  } else {
  	createBlackPlot(res, obsGrid,titleString)
  }
  par(pty=oldpty)
}

createBlackPlot = function(res, obsGrid,titleString ){
#  qpoints = c();
#  rpoints = c();  
#  for(i in 1:length(obsGrid)){
#    idx = which(res[i,] > 0)
#    qpoints = c(qpoints,rep(obsGrid[i], length(idx)))
#    rpoints = c(rpoints,obsGrid[idx])
#  }
#  points(qpoints, rpoints, pch = '.')
  # image(res, col=c('white','black'), axes=FALSE, xlab = 'Observed time grid', ylab = 'Observed time grid', main = titleString)
  
  u1 = as.vector(res)
  u2 = as.vector(t(res))
  t1 = rep(obsGrid, times = length(obsGrid) )
  t2 = rep(obsGrid, each = length(obsGrid)) 
  plot(t1[u1 != 0], t2[u2 !=0] , xlab = 'Observed time grid', ylab = 'Observed time grid', main = titleString, pch = 19, cex =0.33 )
 # axis(1, obsGrid[round(seq(1,length(obsGrid), length.out=11))], obsGrid[round(seq(1,length(obsGrid), length.out=11))],col.axis="black")
 # axis(2, obsGrid[round(seq(1,length(obsGrid), length.out=11))], obsGrid[round(seq(1,length(obsGrid), length.out=11))],col.axis="black")
  
}

createColorPlot = function(res, obsGrid,titleString, ... ){
  res[res > 4] = 4;
  # resVals = sort(unique(as.vector(res)));
  # colPalette = c('white', 'black', 'blue', 'green', 'red')
  # resColPalt = colPalette[resVals+1]
  # image(res, col= resColPalt, axes=FALSE, xlab = 'Observed support points', ylab = 'Observed support points', main = titleString)
  
  notZero <- res != 0
  nnres <- res[notZero]
  ddd <- list(...)
  
  colVec <- c(`1`='black', `2`='blue', `3`='green', `4`='red')
  if (!is.null(ddd[['col']]))
    colVec[] <- ddd[['col']]
    
  pchVec <- rep(19, length(colVec))
  names(pchVec) <- names(colVec)
  if (!is.null(ddd[['pch']]))
    pchVec[] <- ddd[['pch']]
    
  cexVec <- seq(from=0.3, by=0.1, length.out=length(colVec))
  names(cexVec) <- names(colVec)
  if (!is.null(ddd[['cex']]))
    cexVec[] <- ddd[['cex']]
    
  if (!is.null(ddd[['xlab']]))
    xlab <- ddd[['xlab']]
  else
    xlab <- 'Observed time grid'
    
  if (!is.null(ddd[['ylab']]))
    ylab <- ddd[['ylab']]
  else
    ylab <- 'Observed time grid'
  
  t1 = rep(obsGrid, times = length(obsGrid))
  t2 = rep(obsGrid, each = length(obsGrid)) 
  plot(t1[notZero], t2[notZero], col= colVec[nnres], xlab=xlab, ylab=ylab, main = titleString, pch = pchVec[nnres], cex=cexVec[nnres] )
  
  if (!identical(unique(nnres), 1))
    legend('right', colVec, pch = pchVec, col=colVec, pt.cex=cexVec, title = 'Count',bg='white' )
}





