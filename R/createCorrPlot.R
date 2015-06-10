# This function creates the correlation surface plot based on the
# results from FPCA() or FPCder()
######
# Input
######
#  yy : returned object from FPCA().
#  corrPlotType: 'Smoothed': plot the smoothed cov surface 
#                'Fitted': plot the fitted cov surface
#                'Raw': plot the raw covariance scatter plot
#  isInteractive: TRUE: interactive plot
#                 FALSE: printable plot
#  ...: passed into persp3d, persp3D, plot3d or points3D for plotting options

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

  # Check if rgl is installed
  if(isInteractive == TRUE && is.element('rgl', installed.packages()[,1]) == FALSE){
    isInteractive = FALSE
    warning("Interactive plot requires package 'rgl', isInteractive set to be FALSE!")
  }
  
  yname = yy$yname
  workGrid = yy$workGrid
  if(corrPlotType == 'Smoothed'){
    covSurf = yy$smoothedCov # smoothed covariance matrix
  } else if(corrPlotType == 'Fitted'){
    covSurf = yy$fittedCov # fitted covariance matrix
  } else if(corrPlotType == 'Raw'){
    covSurf = yy$rawCov # raw covariance object
  }
  
  if(isInteractive){ # Interactive Plot
    if(corrPlotType == 'Fitted'){
      persp3d(workGrid, workGrid, covSurf, col = 'blue',
        xlab = 't', ylab = 't', zlab = 'Fitted Correlation Surface',
        main = paste('Fitted correlation surface for function', yname))
    } else if(corrPlotType == 'Smoothed'){
      persp3d(workGrid, workGrid, covSurf, col = 'blue',
        xlab = 't', ylab = 't', zlab = 'Smoothed Correlation Surface',
        main = paste('Smoothed correlation surface for function', yname))
    } else if(corrPlotType == 'Raw'){
      tPairs = covSurf$tPairs
      cxxn = covSurf$cxxn
      plot3d(tPairs[,1], tPairs[,2], cxxn, pch = 19, col = 'blue',
        xlab = 't', ylab = 't', zlab = 'Raw Covariance Scatter Plot',
        main = paste('Raw covariance scatter plot for function', yname))
    }
  } else { # Static plot
    if(corrPlotType == 'Fitted'){
      persp3D(workGrid, workGrid, covSurf,
        xlab = 't', ylab = 't', zlab = 'Fitted Correlation Surface',
        main = paste('Fitted correlation surface for function', yname))
    } else if(corrPlotType == 'Smoothed'){
      persp3D(workGrid, workGrid, covSurf,
        xlab = 't', ylab = 't', zlab = 'Smoothed Correlation Surface',
        main = paste('Smoothed correlation surface for function', yname))      
    } else if(corrPlotType == 'Raw'){
      tPairs = covSurf$tPairs
      cxxn = covSurf$cxxn
      points3D(tPairs[,1], tPairs[,2], cxxn, pch = 19,
        xlab = 't', ylab = 't', zlab = 'Raw Covariance Scatter Plot',
        main = paste('Raw covariance scatter plot for function', yname))
    }
  }

}