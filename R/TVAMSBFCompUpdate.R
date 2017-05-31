#####
##### smooth backfitting for a component function
#####

##### input variables: 
#####   f: current SBF estimator of component functions at each estimation points (N*d matrix)
#####   ind: index of updating component during SBF algorithm (scalar)
#####   fNW: marginal regression function kernel estimators at each estimation points (N*d matrix)
#####   Y: response observation points (n-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   X: covariate observation grid (n*d matrix)
#####   h: bandwidths (d-dim. vector)
#####   K: kernel function (function object, default is the Epanechnikov kernel)
#####   supp: supports of estimation interested (d*2 matrix)
#####   MgnJntDensity: evaluated values of marginal and 2-dim. joint densities (2-dim. list, referred to the output of 'MgnJntDensity')

##### output:
#####   updated smooth backfitting component functions for a designated component (N*d matrix)

TVAMSBFCompUpdate<-function(f,ind,fLL,Y,T,X,t,x,h0=NULL,h=NULL,K='epan',supp0=NULL,supp=NULL,MgnJntDens){
  
  N <- nrow(x)
  d <- ncol(x)
  n <- length(Y)
  M <- length(t)
  
  if (K!='epan') {
    message('Epanechnikov kernel is only supported currently. It uses Epanechnikov kernel automatically')
    K<-'epan'
  }
  if (is.null(supp0)==TRUE) {
    supp0 <- c(0,1)
  }
  if (is.null(supp)==TRUE) {
    supp <- matrix(rep(c(0,1),d),ncol=2,byrow=TRUE)
  }
  if (is.null(h0)==TRUE) {
    h0 <- 0.25*n^(-1/5)*(supp0[2]-supp0[1])
  }
  if (is.null(h)==TRUE) {
    h <- rep(0.25*n^(-1/5),d)*(supp[,2]-supp[,1])
  }
  
  tmp <- fLL$f0
  f0LL <-array(0,c(M,N,3))
  for (l in 1:N) {
    f0LL[,l,] <- tmp
  }
  
  fjLL <- fLL$f
  
  j <- ind
  
  tmp1 <- tmp2 <- array(0,c(M,N,3))
  if (j==1) {
    
    for (k in (j+1):d) {
      tmp2 <- tmp2+TVAMCondProjection(f,c(k,j),t,x,T,X,MgnJntDens)
    }  
    
    f[,,,j] <- fjLL[,,,j]-f0LL-tmp2
    
  } 
  if (j>1 && j<d) {
    
    for (k in 1:(j-1)) {
      tmp1 <- tmp1+TVAMCondProjection(f,c(k,j),t,x,T,X,MgnJntDens)
    }
    
    for (k in (j+1):d) {
      tmp2 <- tmp2+TVAMCondProjection(f,c(k,j),t,x,T,X,MgnJntDens)
    } 
    
    f[,,,j] <- fjLL[,,,j]-f0LL-tmp1-tmp2
    
  } 
  if (j==d) {
    
    for (k in 1:(d-1)) {
      tmp1 <- tmp1+TVAMCondProjection(f,c(k,j),t,x,T,X,MgnJntDens)
    }  
    
    f[,,,j] <- fjLL[,,,j]-f0LL-tmp1 
    
  }
  
  f[,,,j] <- TVAMCompFntCent(f,j,t,x,MgnJntDens)

  return(f)
  
}
