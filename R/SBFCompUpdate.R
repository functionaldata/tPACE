#####
##### smooth backfitting for a component function
#####

##### input variables: 
#####   f: current SBF estimator of component functions at each estimation points (N*d matrix)
#####   ind_up: index of updating component during SBF algorithm (scalar)
#####   Y: response observation points (n-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   X: covariate observation grid (n*d matrix)
#####   h: bandwidths (d-dim. vector)
#####   K: kernel function (function object, default is the Epanechnikov kernel)
#####   supp: supports of estimation interested (d*2 matrix)
#####   MgnJntDensity: evaluated values of marginal and 2-dim. joint densities (2-dim. list, referred to the output of 'MgnJntDensity')

##### output:
#####   updated smooth backfitting component functions for a designated component (N*d matrix)

SBFCompUpdate<-function(f,ind_up,f_nw,Y,X,x,h=NULL,K=NULL,supp=NULL,MgnJntDens){
  
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
  
  #f_nw<-NWMgnReg(Y,x,X,h,Kh,supp)
  
  tmp_index<-rep(1,n)
  for(l in 1:d){
    tmp_index<-tmp_index*dunif(X[,l],supp[l,1],supp[l,2])*(supp[l,2]-supp[l,1])
  }
  tmp_index<-which(tmp_index==1)
  
  Y_mean<-sum(Y[tmp_index])/length(Y)/p_0(X)    
  
  j<-ind_up
  
  tmp1<-tmp2<-0
  if(j==1){
    for(k in (j+1):d){
      tmp2<-tmp2+CondProjection(f,c(k,j),x,X,MgnJntDens)
    }  
    
    f[,j]<-f_nw[,j]-Y_mean-tmp2 
  }
  
  if(j>1 && j<d){
    for(k in 1:(j-1)){
      tmp1<-tmp1+CondProjection(f,c(k,j),x,X,MgnJntDens)
    }
    for(k in (j+1):d){
      tmp2<-tmp2+CondProjection(f,c(k,j),x,X,MgnJntDens)
    }          
    f[,j]<-f_nw[,j]-Y_mean-tmp1-tmp2
  }
  
  if(j==d){
    for(k in 1:(d-1)){
      tmp1<-tmp1+CondProjection(f,c(k,j),x,X,MgnJntDens)
    }    
    f[,d]<-f_nw[,d]-Y_mean-tmp1 
  }
  
  f[,j]<-CompFntCent(f,j,x,MgnJntDens)

  return(f)
  
}


#asdf<-MgnJntDensity(x, X, h)
#f<-f_nw<-qwer

#f0<-f
#for(j in 1:4){
#  f<-SBFCompUpdate(f,j,Y,X,x,h,MgnJntDens=asdf)
#}

#attach(mtcars)
#par_tmp<-par(no.readonly=T)
#par(mfrow=c(2,2))
#for(j in 1:4){
#  plot(x[,j],f[,j],type='l',ylim=c(-2,2))
#  points(x[,j],f_nw[,j],type='l',col=2)
#}
#par(par_tmp)
#detach(mtcars)
