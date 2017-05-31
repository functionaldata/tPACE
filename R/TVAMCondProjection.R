#####
##### conditional projection
#####

##### input variables: 
#####   f: evaluated values of component functions at estimation grid (M*N*3*d matrix)
#####   kj: index of kernel estimation for 2-dim. joint density (2-dim. vector)
#####   t: estimation grid (M-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   T: covariate observation grid (n-dim. list M-dim. vectors)
#####   X: covariate observation grid (n*d matrix)
#####   MgnJntDensity: evaluated values of marginal and 2-dim. joint densities (3-dim. list, referred to the output of 'MgnJntDensity')

##### output:
#####   conditional projection of the k-th component function on the j-th component function space (M*N*3 array)


TVAMCondProjection <- function(f, kj, t, x, T, X, MgnJntDens){
  
  N <- nrow(x)
  d <- ncol(x)
  n <- length(X)
  M <- length(t)
  
  k <- kj[1]
  j <- kj[2]
  
  xj <- x[,j]
  xk <- c()
  
  fk <- f[,,,k]
  if (ncol(fk[,,1])==n) {
    xk <- X[,k]
  } else {
    xk <- x[,k]
  }
  
  qj <- MgnJntDens$qArrMgn[,,,,j]
  qkj <- MgnJntDens$qArrJnt[,,,,,k,j]
  
  fTmp <- array(0,c(M,N,3))
  for (m in 1:M) {
    for (lj in 1:N) {
      
      tmp1 <- matrix(0,nrow=N,ncol=3)
      for (lk in 1:N) {
        tmp1[lk,] <- t(qkj[m,lk,lj,,])%*%fk[m,lk,]
      }
      
      tmp2 <- c()
      tmp2[1] <- trapzRcpp(sort(xk),tmp1[order(xk),1])
      tmp2[2] <- trapzRcpp(sort(xk),tmp1[order(xk),2])
      tmp2[3] <- trapzRcpp(sort(xk),tmp1[order(xk),3])
      
      fTmp[m,lj,] <- solve(qj[m,lj,,])%*%tmp2
      
    }
  }
  
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
# fLL<-TVAMLLMgnReg(Y, t, x, T, X, h0=h0, h=h)
# f<-fLL$f
# 
# f1<-f[,,,1]
# f2<-f[,,,2]
# # f1<-CompFntCent(fLL,1,t,x,MgnJntDens)
# # f2<-CompFntCent(fLL,2,t,x,MgnJntDens)
# #
# # f<-array(cbind(f1,f2),c(31,31,3,2))
# 
# MgnJntDens<-TVAMMgnJntDensity(t, x, T, X, h0=h0, h=h)
# 
# f1Cond2<-CondProjection(f, c(1,2), t, x, T, X, MgnJntDens)
# f2Cond1<-CondProjection(f, c(2,1), t, x, T, X, MgnJntDens)
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
# # example of centering
# par(mfrow=c(2,2))
# persp(t,x0,f1[,,1],theta=-30,phi=30, col = surf.colors(f1[,,1]), lty=2, border=NA, shade=0.1, ticktype='detailed',zlim=c(-1.5,1.5))
# persp(t,x0,f2[,,1],theta=-30,phi=30, col = surf.colors(f2[,,1]), lty=2, border=NA, shade=0.1, ticktype='detailed',zlim=c(-1.5,1.5))
# persp(t,x0,f1Cond2[,,1],theta=-30,phi=30, col = surf.colors(f1Cond2[,,1]), lty=2, border=NA, shade=0.1, ticktype='detailed',zlim=c(-1.5,1.5))
# persp(t,x0,f2Cond1[,,1],theta=-30,phi=30, col = surf.colors(f2Cond1[,,1]), lty=2, border=NA, shade=0.1, ticktype='detailed',zlim=c(-1.5,1.5))
