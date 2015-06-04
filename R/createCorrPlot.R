# This function creates the correlation surface plot based on the
# results from FPCA() or FPCder()
######
# Input
######
#  yy : returned object from FPCA().

createCorrPlot = function(yy){
  no_opt = length(yy$eigVal)
  if(no_opt == NULL){
    warning('Warning: Input is not a valid FPCA or FPCder output.')
  } else if(no_opt == 1){
    warning('Warning: Correlation surface is not available when only one principal component is used.')
  } else {
    yname = yy$yname
    regGrid = yy$regGrid
    xcorr = yy$xcorr
    
    plot3D(regGrid, regGrid, xcorr,
      xlab = 't', ylab = 't', zlab = 'Correlation Surface',
      main = paste('Fitted correlation surface for function', yname))
  }
}