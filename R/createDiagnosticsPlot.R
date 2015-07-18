
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

