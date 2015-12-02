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

createDesignPlot = function(t, obsGrid = NULL, isColorPlot=TRUE, noDiagonal=TRUE, ...){
  
  if( class(t) != 'list'){
    stop("You do need to pass a list argument to 'createDesignPlot'!");
  }
  if( is.null(obsGrid)){
    obsGrid = sort(unique(unlist(t)))
  }
  
  args1 <- list( main="Design Title", xlab= 'Observed time grid', ylab= 'Observed time grid', col='black')
  inargs <- list(...)
  args1[names(inargs)] <- inargs 
  
  res = designPlotCount(t, obsGrid, noDiagonal, isColorPlot)
    
  oldpty <- par()[['pty']]
  par(pty="s")
  if(isColorPlot == TRUE){
    createColorPlot(res, obsGrid, args1)
  } else {
    createBlackPlot(res, obsGrid, args1)
  }
  par(pty=oldpty)
  
}

createBlackPlot = function(res, obsGrid, args1){
 
  u1 = as.vector(res)
  u2 = as.vector(t(res))
  t1 = rep(obsGrid, times = length(obsGrid) )
  t2 = rep(obsGrid, each = length(obsGrid)) 
  do.call( plot, c(args1, list(  pch = 19, cex =0.33, x = t1[u1 != 0], y = t2[u2 !=0] ) ) )  
  
}

createColorPlot = function(res, obsGrid, args1){
  
  res[res > 4] = 4;
  notZero <- res != 0
  nnres <- res[notZero]

  # The following three are too involved for a user to define outside 
  # the function as they would need to read the function before-hand,
  # no need to allow control over them from outside
  
  args1$col = NULL
  colVec <- c(`1`='black', `2`='blue', `3`='green', `4`='red')
  # if (!is.null(ddd[['col']])){ 
  #  colVec[] <- ddd[['col']]
  # }
  
  pchVec <- rep(19, length(colVec))
  names(pchVec) <- names(colVec)
  # if (!is.null(ddd[['pch']])){
  #   pchVec[] <- ddd[['pch']]
  # }
  
  cexVec <- seq(from=0.3, by=0.1, length.out=length(colVec))
  names(cexVec) <- names(colVec)
  # if (!is.null(ddd[['cex']])){
  #   cexVec[] <- ddd[['cex']]
  # }
   
  t1 = rep(obsGrid, times = length(obsGrid))
  t2 = rep(obsGrid, each = length(obsGrid)) 
  do.call( plot, c(args1, list( x = t1[notZero], y = t2[notZero], col= colVec[nnres],  pch = pchVec[nnres], cex = cexVec[nnres] ) ))
  
  if (!identical(unique(nnres), 1)){
    legend('right', c('1','2','3','4+'), pch = pchVec, col=colVec, pt.cex=1.5, title = 'Count',bg='white' )
  }
  
}


