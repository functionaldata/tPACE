# This function creates the design plot of the data
######
# Input:
###### 
# t:  input time cell array
# obsGrid:  vector of sorted observed time points
# isColorPlot: TRUE: create color plot with color indicating counts
#             FALSE: create black and white plot with dots indicating observed time pairs
# noDiagonal:  TRUE:  remove diagonal time pairs 
#             FALSE:  do not remove diagonal time pairs
# yname:  the name of the variable containing functional observations

# function createDesignPlot(datafile, isColorPlot, noDiagonal, yname)


createDesignPlot = function(t, obsGrid, isColorPlot, noDiagonal, yname){
  res = designPlotCount(t, obsGrid, noDiagonal, isColorPlot)
  # binrawcov is the function in place of designPlotCount
  # or getCount

  # construct plot structure
  plot(obsGrid[1], obsGrid[1], xlim = range(obsGrid) + isColorPlot * 0.25 * c(0,diff(range(obsGrid))),
     ylim = range(obsGrid),
     xlab = 'Observed time grid', ylab = 'Observed time grid',
     main = paste('Design Plot of', yname), col = 'white')

  if(isColorPlot == TRUE){
  	createColorPlot(res, obsGrid)
  } else {
  	createBlackPlot(res, obsGrid)
  }

}

createBlackPlot = function(res, obsGrid){
  for(i in 1:length(obsGrid)){
    idx = which(res[i,] > 0)
    points(rep(obsGrid[i], length(idx)), obsGrid[idx], pch = 19)
  }
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
  legend('right', c('1', '2', '3~5', '>=6'), pch = 19, 
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