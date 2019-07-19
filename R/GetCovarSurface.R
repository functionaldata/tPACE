GetCovarSurface=function(Ly, Lt, optns = list()){
  # Set the options structure members that are still NULL
  optns = SetOptions(Ly, Lt, optns);
  
  scsObj <- GetScsObj(Ly, Lt, optns)
  covSurf <- scsObj$smoothCov
  if(optns$plot){
    if ((!'plot3D' %in% installed.packages()[, ])) {
      stop("CreateCovPlot requires package 'plot3D'")
    }
    workGrid = scsObj$outGrid
    args = list (x = workGrid, y = workGrid, z = covSurf)
    do.call(plot3D::persp3D, c(args, 'Fitted covariance surface')) 
  }
  return(covSurf)
}
