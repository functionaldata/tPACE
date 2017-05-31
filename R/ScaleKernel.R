#####
##### scaling of compactly supported kernel
#####

##### input variables
#####   x: estimation points (N-dim. vector)
#####   X: covariate observation points (n-dim. vector)
#####   h: bandwidth (scalar)
#####   K: kernel function (function object, default is the Epanechnikov kernel)
#####   supp: support of stimation interested (2-dim. vector, default is a closed interval [0,1])

##### output variable: 
#####   evaluated values of scaled kernel at estimation points near observation points (N by n matrix)


### scaled kernel
ScaleKernel <- function(x, X, h=NULL, K='epan',supp=NULL){
  
  N <- length(x)
  n <- length(X)
  
  if (K!='epan') {
    message('Epanechnikov kernel is only supported currently. It uses Epanechnikov kernel automatically')
  }
  if (is.null(supp)==TRUE) {
    supp <- c(0,1)
  }
  if (is.null(h)==TRUE) {
    h <- 0.25*n^(-1/5)*(supp[2]-supp[1])
  }
  
  xTmp <- matrix(rep(x,n),nrow=N)
  XTmp <- matrix(rep(X,N),ncol=n,byrow=TRUE)
  
  Tmp <- xTmp-XTmp
  
  KhTmp <- (3/4)*(1-(Tmp/h)^2)*dunif(Tmp/h,-1,1)*2/h
  
  return(KhTmp)
}