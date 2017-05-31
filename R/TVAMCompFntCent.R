#####
##### centering conponent function by marginal mean
#####

##### input variables: 
#####   f: evaluated values of component functions at estimation grid (M*N*3*d matrix)
#####   j: index of centering for the j-th component function (scalar)
#####   t: estimation grid (M-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   MgnJntDensity: evaluated values of marginal and 2-dim. joint densities (3-dim. list, referred to the output of 'MgnJntDensity')

##### output:
#####   Locally linear marginal regression function kernel estimators at estimation grid (M*N*3 array)


# centering
TVAMCompFntCent <- function(f,j,t,x,MgnJntDens){
  
  N <- nrow(x)
  d <- ncol(x)
  M <- length(t)
  
  fj0 <- f[,,1,j]
  fj1 <- f[,,2,j]
  fj2 <- f[,,3,j]
  
  xj <- x[,j]
  
  q0 <- MgnJntDens$q0
  qArrMgn <- MgnJntDens$qArrMgn
  
  qj0 <- qArrMgn[,,1,1,j]
  qj1 <- qArrMgn[,,2,1,j]
  qj2 <- qArrMgn[,,3,1,j]
  
  # tmp0 <- tmp1 <- tmp2 <- matrix(0,nrow=M,ncol=N)
  # for(m in 1:M) {
  #   tmp0[m,] <- fj0[m,]-trapzRcpp(sort(xj),(fj0[m,]*qj0[m,]/q0[m])[order(xj)])
  #   tmp1[m,] <- fj1[m,]-trapzRcpp(sort(xj),(fj1[m,]*qj1[m,]/q0[m])[order(xj)])
  #   tmp2[m,] <- fj2[m,]-trapzRcpp(sort(xj),(fj2[m,]*qj2[m,]/q0[m])[order(xj)])
  # }
  # 
  # fTmp <- array(cbind(tmp0,tmp1,tmp2),c(M,N,3))
  
  tmp <- matrix(0,nrow=M,ncol=N)
  for(m in 1:M) {
    #tmp0 <- c(fTmp[m,,1]*qj0[m,]/q0[m])
    tmp1 <- c(fj1[m,]*qj1[m,]/q0[m])
    tmp2 <- c(fj2[m,]*qj2[m,]/q0[m])

    tmp[m,] <- fj0[m,]-trapzRcpp(sort(xj),(tmp1+tmp2)[order(xj)])
  }

  fTmp <- array(cbind(tmp,fj1,fj2),c(M,N,3))
  
  return(fTmp)
  
}


# # test example
# M<-41
# N<-41
# n<-200
# d<-2
# 
# t<-seq(0,1,length.out=M)
# x0<-seq(0,1,length.out=N)
# x<-cbind(x0,x0)
# 
# T<-Y<-list()
# X<-matrix(nrow=n,ncol=d)
# g1 <- function(x1) 2*(x1-0.5)
# g2 <- function(x2) sin(2*pi*x2)
# 
# g <- function(u,x) sin(2*pi*u)*g1(x[1]) + cos(2*pi*u)*g2(x[2])
# for(i in 1:n){
#   X[i,]<-c(pnorm(matrix(c(1,0.5,0.5,1),nrow=2)%*%rnorm(2)))#runif(2)
#   T[[i]]<-t
#   Y[[i]]<-g(T[[i]],X[i,])+rnorm(M,0,0.5)
# }
# 
# K<-'epan'
# h0<-0.15
# h<-c(0.15,0.15)
# 
# fLLMgn<-TVAMLLMgnReg(Y, t, x, T, X, h0=h0, h=h)
# fLL<-fLLMgn$f
# 
# MgnJntDens<-TVAMMgnJntDensity(t, x, T, X, h0=h0, h=h)
# 
# f1Cent<-TVAMCompFntCent(fLL,1,t,x,MgnJntDens)[,,1]
# f2Cent<-TVAMCompFntCent(fLL,2,t,x,MgnJntDens)[,,1]
# 
# surf.colors <- function(x, col = heat.colors(200)) {
# 
#   # First we drop the 'borders' and average the facet corners
#   # we need (nx - 1)(ny - 1) facet colours!
#   x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
#               x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
# 
#   # Now we construct the actual colours matrix
#   colors = col[cut(x.avg, breaks = length(col), include.lowest = F)]
# 
#   return(colors)
# }
# 
# f1<-fLL[,,1,1]
# f2<-fLL[,,1,2]
# 
# # example of centering
# par(mfrow=c(2,2))
# persp(t,x0,f1,theta=-30,phi=30, col = surf.colors(f1), lty=2, border=NA, shade=0.1,ticktype='detailed')
# persp(t,x0,f2,theta=-30,phi=30, col = surf.colors(f2), lty=2, border=NA, shade=0.1,ticktype='detailed')
# persp(t,x0,f1Cent,theta=-30,phi=30, col = surf.colors(f1Cent), lty=2, border=NA, shade=0.1,ticktype='detailed')
# persp(t,x0,f2Cent,theta=-30,phi=30, col = surf.colors(f2Cent), lty=2, border=NA, shade=0.1,ticktype='detailed')
# 
# trapzRcpp(sort(x[,1]),f1Cent[1,order(x[,1])])-trapzRcpp(sort(x[,1]),fLL[1,order(x[,1]),2,1])-trapzRcpp(sort(x[,1]),fLL[1,order(x[,1]),3,1])

