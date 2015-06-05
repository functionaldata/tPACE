# This function creates the correlation surface plot based on the
# results from FPCA() or FPCder()
######
# Input
######
#  yy : returned object from FPCA().
#  corrPlotType: 'Smoothed': plot the smoothed cov surface 
#                'Fitted': plot the fitted cov surface
#  isInteractive: TRUE: interactive plot
#                 FALSE: printable plot
#  ...: passed into persp3d or persp3D for plotting options

createCorrPlot = function(yy, corrPlotType = 'Fitted', isInteractive = FALSE){
  # Check if plotting covariance surface for fitted covariance surface is proper
  if(corrPlotType == 'Fitted'){
    no_opt = length(yy$eigVal)
    if(no_opt == NULL){
      warning('Warning: Input is not a valid FPCA or FPCder output.')
      return()
    } else if(no_opt == 1){
      warning('Warning: Correlation surface is not available when only one principal component is used.')
      return()
    }
  }
  
  yname = yy$yname
  regGrid = yy$regGrid
  truncGrid = yy$truncatedRegGrid
  smoothCov = yy$smoothedCov
  fitCov = yy$fittedCov
  
  if(isInteractive){ # Interactive Plot
    if(corrPlotType == 'Fitted'){
      persp3d(truncGrid, truncGrid, fitCov, col = 'blue',
        xlab = 't', ylab = 't', zlab = 'Fitted Correlation Surface',
        main = paste('Fitted correlation surface for function', yname))
    } else if(corrPlotType == 'Smoothed'){
      persp3d(regGrid, regGrid, smoothCov, col = 'blue',
        xlab = 't', ylab = 't', zlab = 'Smoothed Correlation Surface',
        main = paste('Smoothed correlation surface for function', yname))

    }
  } else { # Static plot
    if(corrPlotType == 'Fitted'){
      persp3D(truncGrid, truncGrid, fitCov,
        xlab = 't', ylab = 't', zlab = 'Fitted Correlation Surface',
        main = paste('Fitted correlation surface for function', yname))
    } else if(corrPlotType == 'Smoothed'){
      persp3D(regGrid, regGrid, smoothCov,
        xlab = 't', ylab = 't', zlab = 'Smoothed Correlation Surface',
        main = paste('Smoothed correlation surface for function', yname))      
    }
  }

}