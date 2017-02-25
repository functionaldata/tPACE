#####
##### marginal and 2-dim. joint kernel densities estimators
#####

##### input variables:
#####   ind_j: index of kernel estimation for marginal density (scalar)
#####   ind_kj: index of kernel estimation for 2-dim. joint density (2-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   X: covariate observation grid (n*d matrix)
#####   h: bandwidths (d-dim. vector)
#####   K: kernel function (function object, default is the Epanechnikov kernel)
#####   supp: supports of estimation interested (d*2 matrix)

##### output:
#####   margianl densities at estimation points near observation points (N*d matrix)
#####   2-dim. joint densities at estimation grid near observation grid (N*N*d*d array)


### propertion of non-truncated observation
p_0<-function(X, supp=NULL){ 
  
  n<-nrow(X)
  d<-ncol(X)
  if(is.null(supp)==T){
    supp<-matrix(rep(c(0,1),d),ncol=2,byrow=T)
  }
  
  tmp<-rep(1,n)
  for(j in 1:d){
    tmp<-tmp*dunif(X[,j],supp[j,1],supp[j,2])*(supp[j,2]-supp[j,1])
  }
  
  return(mean(tmp))
}

# marginal density estimation
p_j<-function(ind_j, x, X, h=NULL, K=NULL, supp=NULL){
  
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
  #if(is.null(Kh)==T){
  #  Kh<-Kh
  #}
  
  j<-ind_j

  tmp_index<-rep(1,n)
  for(l in 1:d){
    tmp_index<-tmp_index*dunif(X[,l],supp[l,1],supp[l,2])*(supp[l,2]-supp[l,1])
  }
  index<-which(tmp_index==1)
  p_hat<-apply(NormKernel(x[,j],X[,j],h[j],K,c(supp[j,1],supp[j,2]))[,index],1,'sum')/n

  p_hat<-p_hat/trapzRcpp(sort(x[,j]),p_hat[order(x[,j])])

  return(p_hat/p_0(X,supp))   
}

# 2-dimensional joint density estimation
p_kj<-function(ind_kj, x, X, h=NULL, K=NULL, supp=NULL){
  
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
  #if(is.null(Kh)==T){
  #  Kh<-Kh
  #}
  
  k<-ind_kj[1]
  p_hat_k<-NormKernel(x[,k],X[,k],h[k],K,c(supp[k,1],supp[k,2]))
  
  j<-ind_kj[2]
  p_hat_j<-NormKernel(x[,j],X[,j],h[j],K,c(supp[j,1],supp[j,2]))
  
  tmp_index<-rep(1,n)
  for(l in 1:d){
    tmp_index<-tmp_index*dunif(X[,l],supp[l,1],supp[l,2])*(supp[l,2]-supp[l,1])
  }
  index<-which(tmp_index==1)
  p_hat<-p_hat_k[,index]%*%t(p_hat_j[,index])/n
  
  p_hat<-p_hat/trapzRcpp(sort(x[,j]),p_j(j,x,X,h)[order(x[,j])])/trapzRcpp(sort(x[,k]),p_j(k,x,X,h)[order(x[,k])])
  
  return(p_hat/p_0(X,supp))     
}

# construction of evaluation matrices for marginal and joint densities estimators
MgnJntDensity<-function(x, X, h=NULL, K=NULL, supp=NULL){
  
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
  #if(is.null(Kh)==T){
  #  Kh<-Kh
  #}
  
  p_mat_mgn<-matrix(0,nrow=N,ncol=d)
  p_arr_jnt<-array(0,dim=c(N,N,d,d))
  #cat(paste('Computing all pairs of 1-/2-dim.l marginal/joint density estimators...','\n',sep=''))
  for(j in 1:d){
    #cat(paste('   ',round(j/d,3),'\n',sep=''))
    #cat('\n')
    p_mat_mgn[,j]<-p_j(j,x,X,h,K,supp)
    
    for(k in j:d){ 
      #print(k)
      if(k==j){
        p_arr_jnt[,,k,j]<-diag(p_j(j,x,X,h,K,supp))
      }else{
        p_arr_jnt[,,k,j]<-p_kj(c(k,j),x,X,h,K,supp)
        p_arr_jnt[,,j,k]<-t(p_arr_jnt[,,k,j])
      }
    }
  }
  
  return(list(p_arr_jnt=p_arr_jnt, p_mat_mgn=p_mat_mgn))
}

