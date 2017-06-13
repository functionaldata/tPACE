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

SBFCompUpdate<-function(f,ind,fNW,Y,X,x,h=NULL,K='epan',supp=NULL,MgnJntDens){
  
  N<-nrow(x)
  d<-ncol(x)
  n<-nrow(X)
  
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
  
  #f_nw<-NWMgnReg(Y,x,X,h,Kh,supp)
  
  #tmp <- X >= matrix(supp[, 1], nrow=n, byrow=TRUE) & X <= matrix(supp[, 2], nrow=n, byrow=TRUE)
  # tmp_ind <- apply(tmp, 1, all)
  
  tmpIndex <- rep(1,n)
  for (l in 1:d) {
    tmpIndex <- tmpIndex*dunif(X[,l],supp[l,1],supp[l,2])*(supp[l,2]-supp[l,1])
  }
  tmpIndex <- which(tmpIndex==1)
  
  yMean <- sum(Y[tmpIndex])/length(Y)/P0(X,supp)    
  
  j <- ind
  
  tmp1<-tmp2<-0
  if (j==1) {
    for (k in (j+1):d) {
      tmp2 <- tmp2+CondProjection(f,c(k,j),x,X,MgnJntDens)
    }  
    
    f[,j]<-fNW[,j]-yMean-tmp2 
  } else if (j>1 && j<d) {
    for (k in 1:(j-1)) {
      tmp1 <- tmp1+CondProjection(f,c(k,j),x,X,MgnJntDens)
    }
    for (k in (j+1):d) {
      tmp2 <- tmp2+CondProjection(f,c(k,j),x,X,MgnJntDens)
    }          
    f[,j] <- fNW[,j]-yMean-tmp1-tmp2
  } else if (j==d) {
    for (k in 1:(d-1)) {
      tmp1 <- tmp1+CondProjection(f,c(k,j),x,X,MgnJntDens)
    }    
    f[,d] <- fNW[,d]-yMean-tmp1 
  }
  
  f[,j] <- CompFntCent(f,j,x,MgnJntDens)

  return(f)
  
}
