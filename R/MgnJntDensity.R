#####
##### marginal and 2-dim. joint kernel densities estimators
#####

##### input variables:
#####   j: index of kernel estimation for marginal density (scalar)
#####   kj: index of kernel estimation for 2-dim. joint density (2-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   X: covariate observation grid (n*d matrix)
#####   h: bandwidths (d-dim. vector)
#####   K: kernel function (function object, default is the Epanechnikov kernel)
#####   supp: supports of estimation interested (d*2 matrix)

##### output:
#####   margianl densities at estimation points near observation points (N*d matrix)
#####   2-dim. joint densities at estimation grid near observation grid (N*N*d*d array)


### propertion of non-truncated observation
P0 <- function(X, supp=NULL){ 
  
  n <- nrow(X)
  d <- ncol(X)
  
  if (is.null(supp)==TRUE) {
    supp <- matrix(rep(c(0,1),d),ncol=2,byrow=TRUE)
  }
  
  tmp <- rep(1,n)
  for(j in 1:d){
    tmp <- tmp*dunif(X[,j],supp[j,1],supp[j,2])*(supp[j,2]-supp[j,1])
  }
  
  return(mean(tmp))
}

# marginal density estimation
Pj <- function(j, x, X, h=NULL, K='epan', supp=NULL){
  
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

  tmpIndex <- rep(1,n)
  for(l in 1:d){
    tmpIndex <- tmpIndex*dunif(X[,l],supp[l,1],supp[l,2])*(supp[l,2]-supp[l,1])
  }
  index <- which(tmpIndex==1)
  pHat <- apply(NormKernel(x[,j],X[,j],h[j],K,c(supp[j,1],supp[j,2]))[,index],1,'sum')/n

  pHat <- pHat/trapzRcpp(sort(x[,j]),pHat[order(x[,j])])

  return(pHat/P0(X,supp))   
}

# 2-dimensional joint density estimation
Pkj <- function(kj, x, X, h=NULL, K='epan', supp=NULL){
  
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
  
  k <- kj[1]
  pHatk <- NormKernel(x[,k],X[,k],h[k],K,c(supp[k,1],supp[k,2]))
  
  j <- kj[2]
  pHatj <- NormKernel(x[,j],X[,j],h[j],K,c(supp[j,1],supp[j,2]))
  
  tmpIndex <- rep(1,n)
  for(l in 1:d){
    tmpIndex <- tmpIndex*dunif(X[,l],supp[l,1],supp[l,2])*(supp[l,2]-supp[l,1])
  }
  index <- which(tmpIndex==1)
  pHat <- pHatk[,index]%*%t(pHatj[,index])/n
  
  pHat <- pHat/trapzRcpp(sort(x[,j]),Pj(j,x,X,h,K,supp)[order(x[,j])])/trapzRcpp(sort(x[,k]),Pj(k,x,X,h,K,supp)[order(x[,k])])
  
  return(pHat/P0(X,supp))     
}

# construction of evaluation matrices for marginal and joint densities estimators
MgnJntDensity <- function(x, X, h=NULL, K='epan', supp=NULL){
  
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
  
  pMatMgn <- matrix(0,nrow=N,ncol=d)
  pArrJnt <- array(0,dim=c(N,N,d,d))
  #cat(paste('Computing all pairs of 1-/2-dim.l marginal/joint density estimators...','\n',sep=''))
  for (j in 1:d) {
    #cat(paste('   ',round(j/d,3),'\n',sep=''))
    #cat('\n')
    pMatMgn[,j] <- Pj(j,x,X,h,K,supp)
    
    for (k in j:d) { 
      #print(k)
      if (k==j) {
        #pArrJnt[,,k,j] <- diag(Pj(j,x,X,h,K,supp))
        pArrJnt[,,k,j] <- diag(pMatMgn[,j])
      } else {
        pArrJnt[,,k,j] <- Pkj(c(k,j),x,X,h,K,supp)
        pArrJnt[,,j,k] <- t(pArrJnt[,,k,j])
      }
    }
  }
  
  return(list(pArrJnt=pArrJnt, pMatMgn=pMatMgn))
}

