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
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' createDesignPlot(sampWiener$tList, sort(unique(unlist(sampWiener$tList))))
#' @export

createDesignPlot = function(t, obsGrid = NULL, isColorPlot=FALSE, noDiagonal=TRUE, yname = NULL){
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
  plot(obsGrid[1], obsGrid[1], xlim = range(obsGrid) + isColorPlot * 0.25 * c(0,diff(range(obsGrid))),
     ylim = range(obsGrid),
     xlab = 'Observed time grid', ylab = 'Observed time grid',
     main = titleString, col = 'white')

  if(isColorPlot == TRUE){
  	createColorPlot(res, obsGrid)
  } else {
  	createBlackPlot(res, obsGrid)
  }

}

createBlackPlot = function(res, obsGrid){
  qpoints = c();
  rpoints = c();
  
  for(i in 1:length(obsGrid)){
    idx = which(res[i,] > 0)
    qpoints = c(qpoints,rep(obsGrid[i], length(idx)))
    rpoints = c(rpoints,obsGrid[idx])
  }
  points(qpoints, rpoints, pch = '.')
}

createColorPlot = function(res, obsGrid){
  for(i in 1:length(obsGrid)){
    tmp = res[i,]
    idx = which(tmp > 0)
    if(length(idx) > 0){
      for(j in 1:length(idx)){
        points(obsGrid[i], obsGrid[idx[j]], col = searchCol(tmp[idx[j]]),
          pch = 19)
      }
    }
  }
  legend('right', c('1', '2', '3~5', '>=6'), pch = 14, 
    col = c('red', 'purple', 'green', 'blue'), title = 'Count')
}

searchCol = function(val){
  if(val == 1){
	  col = 'red'
	} else if(val == 2){
	  col = 'purple'
	} else if(val >= 3 && val <= 5){
    col = 'green'
  } else if(val >= 6){
    col = 'blue'
  } else {
    col = 'white'
  }
  return(col)
}
