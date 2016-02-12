#' Functional Principal Component Analysis Diagnostics plot
#' 
#' This function will open a new device if not instructed otherwise.
#'
#' @param ret An FPCA class object returned by FPCA().
#' @param openNewDev A logical specifying if a new device should be opened - default: FALSE
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=FALSE))
#' CreateDiagnosticsPlot(res)
#' @export

CreateDiagnosticsPlot <-function(ret, openNewDev = FALSE){ 
  t = ret$inputData$t
  if(openNewDev){ 
    dev.new(width=6.2, height=6.2, noRStudioGD=TRUE) ; 
  }
  fves = ret$cumFVE
  mu = ret$mu
  obsGrid = ret$obsGrid      
  workGrid = ret$workGrid
 
  par(mfrow=c(2,2))
  ## Make Design plot
  CreateDesignPlot(t)
  
  ## Make Mean trajectory plot
  plot( obsGrid, mu, type='l', xlab='s',ylab='', main='Mean Function', panel.first = grid())   
  
  ## Make Scree plot
  CreateScreePlot(ret);
  
  ## Make Phi plot
  K = ncol(ret$phi);
  k =1;
  if(K>3){
    k = 3;
  } else {
    k = K;
  } # paste(c("First ", as.character(3), " eigenfunctions"),collapse= '')
  matplot(workGrid, ret$phi[,1:k], type='n', main=paste(collapse='', c("First ", as.character(k), " Eigenfunctions"))  , xlab='s', ylab='') 
  grid()
  matlines(workGrid, ret$phi[,1:k] ) 
}

