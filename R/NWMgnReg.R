#####
##### Nadaraya-Watson marginal regression estimation
#####

##### input variables: 
#####   Y: response observation points (n-dim. vector)
#####   index_kj: index of conditional projection for the k-th component function on the j-th component function space (2-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   X: covariate observation grid (n*d matrix)
#####   h: bandwidths (d-dim. vector)
#####   K: kernel function (function object, default is the Epanechnikov kernel)
#####   supp: supports of estimation interested (d*2 matrix)

##### output:
#####   NW marginal regression function kernel estimators at eacj estimation points (N*d matrix)

NWMgnReg<-function(Y, x, X, h=NULL, K=NULL, supp=NULL){
  
  N<-nrow(x)
  d<-ncol(x)
  n<-nrow(X)
  if(is.null(K)==T){
    K<-EpchKer
  }
  if(is.null(supp)==T){
    supp<-matrix(rep(c(0,1),d),ncol=2,byrow=T)
  }
  if(is.null(h)==T){
    h<-rep(0.25*n^(-1/5),d)*(supp[,2]-supp[,1])
  }
  
  f_nw<-matrix(0,nrow=N,ncol=d)
  
  tmp_index<-rep(1,n)
  for(j in 1:d){
    tmp_index<-tmp_index*dunif(X[,j],supp[j,1],supp[j,2])*(supp[j,2]-supp[j,1])
  }
  tmp_index<-which(tmp_index==1)
  
  for(j in 1:d){
    p_hat_j<-NormKernel(x[,j],X[,j],h[j],K,c(supp[j,1],supp[j,2]))
    r_hat_j<-c(p_hat_j[,tmp_index]%*%Y[tmp_index])/length(Y)
    
    p_hat_j<-apply(p_hat_j[,tmp_index],1,'sum')/length(Y)
    
    tmp_ind<-which(p_hat_j!=0)
    
    f_nw[tmp_ind,j]<-r_hat_j[tmp_ind]/p_hat_j[tmp_ind]
  }
  
  return(f_nw)
}



#x<-seq(0,1,length.out=11)
#x<-cbind(x,x,x,x)
#X<-matrix(runif(4*400,-0.01,1.01),ncol=4)
#h<-rep(0.1,4)

#Y<-X[,1]+X[,2]^2+cos(2*pi*X[,3])+cos(4*pi*X[,4])
#qwer<-NWMgnReg(Y, x, X, h)

#attach(mtcars)
#par_tmp<-par(no.readonly=T)
#par(mfrow=c(2,2))
#for(j in 1:4){
# plot(x[,j],qwer[,j],type='l')
#  points(X[,j],Y,cex=0.5)
#}
#par(par_tmp)
#detach(mtcars)

