#####
##### conditional projection
#####

##### input variables: 
#####   f: evaluated values of component functions at estimation grid (N*d matrix)
#####   index_kj: index of conditional projection for the k-th component function on the j-th component function space (2-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   X: covariate observation grid (n*d matrix)
#####   MgnJntDensity: evaluated values of marginal and 2-dim. joint densities (2-dim. list, referred to the output of 'MgnJntDensity')

##### output:
#####   conditional projection of the k-th component function on the j-th component function space (N-dim. vector)


CondProjection<-function(f, ind_kj, x, X, MgnJntDens){
  
  N<-nrow(x)
  n<-nrow(X)
  d<-ncol(X)
  
  k<-ind_kj[1]
  j<-ind_kj[2]
  
  x_j<-x[,j]
  x_k<-c()
  
  f_k<-f[,k]
  if(length(f_k)==n){
    x_k<-X[,k]
  }else{
    x_k<-x[,k]
  }
  
  asdf<-MgnJntDens$p_mat_mgn[,j]
  
  tmp_ind<-which(asdf!=0)
  qwer<-MgnJntDens$p_arr_jnt[,tmp_ind,k,j]
  
  if(length(tmp_ind)>0){
    
    p_hat<-matrix(0,nrow=length(x_k),ncol=length(x_j))
    
    p_hat[,tmp_ind]<-t(t(qwer)/asdf[tmp_ind])
    
    tmp<-c()
    for(l in 1:ncol(p_hat)){
      tmptmp<-f_k*c(p_hat[,l])
      tmp[l]<-trapzRcpp(sort(x_k),tmptmp[order(x_k)])
    }
    
    return(tmp)    
    
  }else{
    return(0)
  }
}
