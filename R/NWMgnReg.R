#####
##### Nadaraya-Watson marginal regression estimation
#####

##### input variables: 
#####   Y: response observation points (n-dim. vector)
#####   kj: index of conditional projection for the k-th component function on the j-th component function space (2-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   X: covariate observation grid (n*d matrix)
#####   h: bandwidths (d-dim. vector)
#####   K: kernel function (function object, default is the Epanechnikov kernel)
#####   supp: supports of estimation interested (d*2 matrix)

##### output:
#####   NW marginal regression function kernel estimators at each estimation point (N*d matrix)

NWMgnReg <- function(Y, x, X, h=NULL, K='epan', supp=NULL){
  
  N <- nrow(x)
  d <- ncol(x)
  n <- nrow(X)
  
  if (K!='epan') {
    message('Epanechnikov kernel is only supported currently. It uses Epanechnikov kernel automatically')
    K<-'epan'
  }
  if (is.null(supp)==TRUE) {
    supp <- matrix(rep(c(0,1),d),ncol=2,byrow=TRUE)
  }
  if (is.null(h)==TRUE) {
    h <- rep(0.25*n^(-1/5),d)*(supp[,2]-supp[,1])
  }
  
  fNW <- matrix(0,nrow=N,ncol=d)
  
  tmpIndex <- rep(1,n)
  for (j in 1:d) {
    tmpIndex <- tmpIndex*dunif(X[,j],supp[j,1],supp[j,2])*(supp[j,2]-supp[j,1])
  }
  tmpIndex <- which(tmpIndex==1)
  
  for (j in 1:d) {
    pHatj <- NormKernel(x[,j],X[,j],h[j],K,c(supp[j,1],supp[j,2]))
    rHatj <- c(pHatj[,tmpIndex]%*%Y[tmpIndex])/length(Y)
    
    pHatj <- apply(pHatj[,tmpIndex],1,'sum')/length(Y)
    
    tmpInd <- which(pHatj!=0)
    
    fNW[tmpInd,j] <- rHatj[tmpInd]/pHatj[tmpInd]
  }
  
  return(fNW)
}

