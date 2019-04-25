#####
##### conditional projection
#####

##### input variables: 
#####   f: evaluated values of component functions at estimation grid (N*d matrix)
#####   kj: index of conditional projection for the k-th component function on the j-th component function space (2-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   X: covariate observation grid (n*d matrix)
#####   MgnJntDensity: evaluated values of marginal and 2-dim. joint densities (2-dim. list, referred to the output of 'MgnJntDensity')

##### output:
#####   conditional projection of the k-th component function on the j-th component function space (N-dim. vector)


CondProjection <- function(f, kj, x, X, MgnJntDens){
  
  N <- nrow(x)
  n <- nrow(X)
  d <- ncol(X)
  
  k <- kj[1]
  j <- kj[2]
  
  xj <- x[,j]
  xk <- c()
  
  fk <- f[,k]
  if (length(fk)==n) {
    xk <- X[,k]
  } else {
    xk <- x[,k]
  }
  
  asdf <- MgnJntDens$pMatMgn[,j]
  
  tmpInd <- which(asdf!=0)
  qwer <- MgnJntDens$pArrJnt[,tmpInd,k,j]
  
  if (length(tmpInd)>0) {
    
    pHat <- matrix(0,nrow=length(xk),ncol=length(xj))
    
    pHat[,tmpInd] <- t(t(qwer)/asdf[tmpInd])
    
    tmp <- c()
    for (l in 1:ncol(pHat)) {
      tmptmp <- fk*c(pHat[,l])
      tmp[l] <- trapzRcpp(sort(xk),tmptmp[order(xk)])
    }
    
    return(tmp)    
    
  } else {
    return(0)
  }
}
