#' Functional Principal Component Analysis Diagnostics plot
#' 
#' FPCA for dense or sparse functional data. 
#' 
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param ret An FPCA class object returned by FPCA().
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createDiagnosticsPlot(sampWiener$tList, res)

createDiagnosticsPlot <-function(t,ret){ 
 dev.new(width=6.2, height=6.2, noRStudioGD=TRUE) ; 
 fves = ret$FVE
 mu = ret$mu
 obsGrid = ret$obsGrid      
 workGrid = ret$workGrid
 
 par(mfrow=c(2,2))
  createDesignPlot(t)
  plot( obsGrid, mu, type='l', xlab='s',ylab='', main='Mean')  
  grid()
  createScreePlot(fves);
  K = length(fves);
  k =1;
  if(K>3){
    k = 3;
  } else {
    k = K;
  }
  matplot(workGrid, ret$phi[,1:k], type='l', main="First Eigenfunctions", xlab='s', ylab='') 
}

