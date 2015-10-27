#' Create the covariance surface plot based on the results from FPCA() or FPCder().
#'
#' This function will open a new device if not instructed otherwise.
#'
#' @param yy returned object from FPCA().
#' @param covPlotType a string specifying the type of covariance surface to be plotted:
#'                     'Smoothed': plot the smoothed cov surface 
#'                     'Fitted': plot the fitted cov surface
#'                     'Raw': plot the raw covariance scatter plot
#' @param isInteractive an option for interactive plot:
#'                      TRUE: interactive plot; FALSE: printable plot
#' @param ... other arguments passed into persp3d, persp3D, plot3d or points3D for plotting options
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createCovPlot(res)
#' @export

createCovPlot = function(yy, covPlotType = 'Fitted', isInteractive = FALSE){
  # Check if plotting covariance surface for fitted covariance surface is proper
  if(covPlotType == 'Fitted'){
    no_opt = length(yy$lambda)
    if(length(no_opt) == 0){
      warning('Warning: Input is not a valid FPCA or FPCder output.')
      return()
    } else if(no_opt == 1){
      warning('Warning: Covariance surface is not available when only one principal component is used.')
      return()
    }
  }

  # Check if rgl is installed
  if(isInteractive == TRUE && is.element('rgl', installed.packages()[,1]) == FALSE){
    isInteractive = FALSE
    warning("Interactive plot requires package 'rgl', isInteractive set to be FALSE!")
  }
  
  if (is.null(yy$yname)){
    yname = NULL
  } else {
    yname = paste(' for function',  yy$yname)
  }
  
  
  workGrid = yy$workGrid
  if(covPlotType == 'Smoothed'){
    covSurf = yy$smoothedCov # smoothed covariance matrix
  } else if(covPlotType == 'Fitted'){
    covSurf = yy$fittedCov # fitted covariance matrix
  } else {
    warning("Covariance plot type no recognised; using default ('Fitted').")
    covSurf = yy$fittedCov 
    covPlotType = 'Fitted'
  }

 # if(covPlotType == 'Raw'){
 #   covSurf = yy$rawCov # raw covariance object
 # }
  
  if(isInteractive){ # Interactive Plot
    if(covPlotType == 'Fitted'){
      rgl::persp3d(workGrid, workGrid, covSurf, col = 'blue',
        xlab = 't', ylab = 't', zlab = 'Fitted Covariance Surface',
        main = paste('Fitted covariance surface', yname))
    } else if(covPlotType == 'Smoothed'){
      rgl::persp3d(workGrid, workGrid, covSurf, col = 'blue',
        xlab = 't', ylab = 't', zlab = 'Smoothed Covariance Surface',
        main = paste('Smoothed covariance surface', yname))
    }
#  else if(covPlotType == 'Raw'){
#      tPairs = covSurf$tPairs
#      cxxn = covSurf$cxxn
#      rgl::plot3d(tPairs[,1], tPairs[,2], cxxn, pch = 19, col = 'blue',
#        xlab = 't', ylab = 't', zlab = 'Raw Covariance Scatter Plot',
#        main = paste('Raw covariance scatter plot', yname))
#    }
  } else { # Static plot
    if(covPlotType == 'Fitted'){
      persp3D(workGrid, workGrid, covSurf,
        xlab = 't', ylab = 't', zlab = 'Fitted Covariance Surface',
        main = paste('Fitted covariance surface', yname))
    } else if(covPlotType == 'Smoothed'){
      persp3D(workGrid, workGrid, covSurf,
        xlab = 't', ylab = 't', zlab = 'Smoothed Covariance Surface',
        main = paste('Smoothed covariance surface', yname))      
    } 
#  else if(covPlotType == 'Raw'){
#      tPairs = covSurf$tPairs
#      cxxn = covSurf$cxxn
#      points3D(tPairs[,1], tPairs[,2], cxxn, pch = 19,
#        xlab = 't', ylab = 't', zlab = 'Raw Covariance Scatter Plot',
#        main = paste('Raw covariance scatter plot', yname))
#    }
  }

}
