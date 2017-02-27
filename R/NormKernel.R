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


### Epanechnikov kernel
EpchKer<-function(t)  (3/4)*(1-t^2)*dunif(t,-1,1)*2

### scaled kernel
Kh<-function(x, X, h=NULL, K=NULL){
  N<-length(x)
  n<-length(X)
  if(is.null(K)==T){
    K<-EpchKer
  }else{
    cat(paste('Epanechnikov kernel is only supported currently.'))
  }
  if(is.null(h)==T){
    h<-0.25*n^(-1/5)*(supp[2]-supp[1])
  }
  
  xTmp<-matrix(rep(x,n),nrow=N)
  XTmp<-matrix(rep(X,N),ncol=n,byrow=T)
  
  return(K((xTmp-XTmp)/h)/h)
}

### normalized kernel
NormKernel<-function(x, X, h=NULL, K=NULL, supp=NULL){
  
  N<-length(x)
  n<-length(X)
  if(is.null(K)==T){
    K<-EpchKer
  }else{
    cat(paste('Epanechnikov kernel is only supported currently.'))
  }
  if(is.null(supp)==T){
    supp<-c(0,1)
  }
  if(is.null(h)==T){
    h<-0.25*n^(-1/5)*(supp[2]-supp[1])
  }
  #if(is.null(Kh)==T){
  #  Kh<-Kh
  #}
  
  numer<-Kh(x,X,h,K=K)
  
  ind1<-which(dunif(X,supp[1],supp[2])==0)
  numer[,ind1]<-0
  
  denom<-c()
  for(i in 1:n){
    denom[i]<-trapzRcpp(sort(x),numer[order(x),i])
  }
  #denom<-apply(numer[x_order,],2,FUN='trapzRcpp',X=sort(x))
  
  ind2<-which(denom==0)
  
  NormKernelTmp<-numer/matrix(rep(denom,N),nrow=N,byrow=T)
  NormKernelTmp[,ind2]<-0
  
  if(min(nrow(NormKernelTmp),ncol(NormKernelTmp))==1){
    return(c(NormKernelTmp))
  }else{
    return(NormKernelTmp)
  }
}

