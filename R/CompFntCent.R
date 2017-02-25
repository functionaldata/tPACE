#####
##### centering conponent function by marginal mean
#####

##### input variables: 
#####   f: evaluated values of component functions at estimation grid (N*d matrix)
#####   ind_j: index of centering for the j-th component function (scalar)
#####   x: estimation grid (N*d matrix)
#####   MgnJntDensity: evaluated values of marginal and 2-dim. joint densities (2-dim. list, referred to the output of 'MgnJntDensity')

##### output:
#####   NW marginal regression function kernel estimators at estimation grid (N*d matrix)


# centering
CompFntCent<-function(f,ind_j,x,MgnJntDens){
  
  j<-ind_j
  
  f_j<-f[,j]
  x_j<-x[,j]
  
  p_mat_mgn<-MgnJntDens$p_mat_mgn
  
  tmp1<-p_mat_mgn[,j]
  tmp<-f_j-trapzRcpp(sort(x_j),(f_j*tmp1)[order(x_j)])
  
  return(tmp)
  
}




#attach(mtcars)
#par_tmp<-par(no.readonly=T)
#par(mfrow=c(2,2))
#for(j in 1:4){
#  plot(x[,j],CompFntCent(qwer,j,x,MgnJntDens=asdf),type='l')
#  points(X[,j],Y,cex=0.5)
#}
#par(par_tmp)
#detach(mtcars)