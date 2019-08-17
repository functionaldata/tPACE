#####
##### normalization of compactly supported kernel (Mammen et al., 1999)
#####

##### input variables
#####   x: estimation points (N-dim. vector)
#####   X: covariate observation points (n-dim. vector)
#####   h: bandwidth (scalar)
#####   K: kernel function (function object, default is the Epanechnikov kernel)
#####   supp: support of stimation interested (2-dim. vector, default is a closed interval [0,1])

##### output variable: 
#####   evaluated values of normalized kernel at estimation points near observation points (N by n matrix)



### normalized kernel
NormKernel <- function(x, X, h=NULL, K='epan', supp=NULL){
  
  N <- length(x)
  n <- length(X)
  
  if (K!='epan') {
    message('Epanechnikov kernel is the default choice')
    K<-'epan'
  }
  if (is.null(supp)==TRUE) {
    supp <- c(0,1)
  }
  if (is.null(h)==TRUE) {
    h <- 0.25*n^(-1/5)*(supp[2]-supp[1])
  }
  
  numer <- ScaleKernel(x,X,h,K=K,supp=supp)
  
  ind1 <- which(dunif(X,supp[1],supp[2])==0)
  numer[,ind1] <- 0
  
  denom <- c()
  for (i in 1:n) {
    denom[i] <- trapzRcpp(sort(x),numer[order(x),i])
  }
  #denom <- apply(numer[x_order,],2,FUN='trapzRcpp',X=sort(x))
  
  ind2 <- which(denom==0)
  
  NormKernelTmp <- numer/matrix(rep(denom,N),nrow=N,byrow=TRUE)
  NormKernelTmp[,ind2] <- 0
  
  if (min(nrow(NormKernelTmp),ncol(NormKernelTmp))==1) {
    return(c(NormKernelTmp))
  } else {
    return(NormKernelTmp)
  }
}

