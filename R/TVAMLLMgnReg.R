#####
##### Locally linear marginal regression estimation for time-varying additive models
#####

##### input variables: 
#####   Y: response observation points (n-dim. list M-dim. vectors)
#####   t: estimation grid (M-dim. vector)
#####   x: estimation grid (N*d matrix)
#####   T: covariate observation grid (n-dim. list M-dim. vectors)
#####   X: covariate observation grid (n*d matrix)
#####   h0: bandwidths (d-dim. vector)
#####   h: bandwidths (d-dim. vector)
#####   K: kernel function (function object, default is the Epanechnikov kernel)
#####   supp0: support of estimation interested (2-dim. vector)
#####   supp: supports of estimation interested (d*2 matrix)

##### output:
#####   f0LL: Locally linear marginal regression function kernel estimators of function values and slopes at each estimation points (M*2 matrix)
#####   fjLL: Locally linear marginal regression function kernel estimators of function values and slopes at each estimation points (M*N*3*d array)


TVAMLLMgnReg <- function(Y, t, x, T, X, h0=NULL, h=NULL, K='epan', supp0=NULL, supp=NULL){
  
  N <- nrow(x)
  d <- ncol(x)
  n <-length(Y)
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
  
  vecY <- vecT <- c()
  tmpRep <- c()
  for (i in 1:n) {
    vecY<-c(vecY,Y[[i]])
    vecT<-c(vecT,T[[i]])
    tmpRep[i] <- length(T[[i]])
  }
  
  # locally linear marginal regression for time
  f0LL <- matrix(0,nrow=M,ncol=3)
  
  v0 <- (matrix(rep(vecT,M), nrow=M, ncol=length(vecT), byrow=TRUE)  - matrix(rep(t,length(vecT)), nrow=M, ncol=length(vecT)))/h0
  
  qHat00 <- NormKernel(t,vecT,h0,K,supp0)
  qHat01 <- qHat00*v0
  qHat02 <- qHat00*v0^2
  
  rHat00 <- qHat00%*%vecY/length(vecY)
  rHat01 <- qHat01%*%vecY/length(vecY)
  
  qHat011 <- apply(qHat00,1,'mean')
  qHat012 <- apply(qHat01,1,'mean')
  qHat022 <- apply(qHat02,1,'mean')
  
  for (m in 1:M) {
    rHat0 <- c(rHat00[m], rHat01[m])
    qHat0 <- matrix(c(qHat011[m],qHat012[m],
                      qHat012[m],qHat022[m]),nrow=2,ncol=2)
  
    f0LL[m,1:2] <- solve(qHat0)%*%rHat0
  }  
  
  # locally linear marginal regression for time and covariate
  fLL <- array(0,dim=c(M,N,3,d))
  
  for (j in 1:d) {
    #vecXj <- rep(X[,j],rep(M,n))
    vecXj <- rep(X[,j],tmpRep)
    
    vj <- (matrix(rep(vecXj,N), nrow=N, ncol=length(vecXj), byrow=TRUE)  - matrix(rep(x[,j],length(vecXj)), nrow=N, ncol=length(vecXj)))/h[j]
  
    normKernelXj <- NormKernel(x[,j],vecXj,h[j],K,supp[j,])
    
    qHatj11 <- qHat00%*%t(normKernelXj)/length(vecY)
    qHatj12 <- (qHat00*v0)%*%t(normKernelXj)/length(vecY)
    qHatj13 <- qHat00%*%t(normKernelXj*vj)/length(vecY)
    qHatj22 <- (qHat00*v0^2)%*%t(normKernelXj)/length(vecY)
    qHatj23 <- (qHat00*v0)%*%t(normKernelXj*vj)/length(vecY)
    qHatj33 <- qHat00%*%t(normKernelXj*vj^2)/length(vecY)

    rHatj0 <- qHat00%*%Matrix::sparseMatrix(1:length(vecY),1:length(vecY),x=vecY)%*%t(normKernelXj)/length(vecY)
    rHatj1 <- (qHat00*v0)%*%Matrix::sparseMatrix(1:length(vecY),1:length(vecY),x=vecY)%*%t(normKernelXj)/length(vecY)
    rHatj2 <- qHat00%*%Matrix::sparseMatrix(1:length(vecY),1:length(vecY),x=vecY)%*%t(normKernelXj*vj)/length(vecY)
    
    # rHatj0 <- qHat00%*%diag(vecY)%*%t(normKernelXj)/length(vecY)
    # rHatj1 <- (qHat00*v0)%*%diag(vecY)%*%t(normKernelXj)/length(vecY)
    # rHatj2 <- qHat00%*%diag(vecY)%*%t(normKernelXj*vj)/length(vecY)
  
    for (m in 1:M) {
      for (k in 1:N) {
        rHatj <- c(rHatj0[m,k], rHatj1[m,k], rHatj2[m,k])
        qHatj <- matrix(c(qHatj11[m,k],qHatj12[m,k],qHatj13[m,k],
                          qHatj12[m,k],qHatj22[m,k],qHatj23[m,k],
                          qHatj13[m,k],qHatj23[m,k],qHatj33[m,k]),nrow=3,ncol=3)
    
        fLL[m,k,,j] <- solve(qHatj)%*%rHatj
      }
    }
  }
  
  return(list(f0=f0LL, f=fLL))
  
}





# 
# # test example
# M<-41
# N<-41
# n<-500
# d<-2
# 
# t<-seq(0,1,length.out=M)
# x0<-seq(0,1,length.out=N)
# x<-cbind(x0,x0)
# 
# T<-Y<-list()
# X<-matrix(nrow=n,ncol=d)
# for(i in 1:n){
#   X[i,]<-runif(2)
#   T[[i]]<-t
#   Y[[i]]<-cos(2*pi*T[[i]])*X[i,1]+rnorm(M,0,0.5)
# }
# 
# K<-'epan'
# h0<-0.1
# h<-c(0.15,0.25)
# 
# surf.colors <- function(x, col = heat.colors(100)) {
# 
#   # First we drop the 'borders' and average the facet corners
#   # we need (nx - 1)(ny - 1) facet colours!
#   x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
#               x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
# 
#   # Now we construct the actual colours matrix
#   colors = col[cut(x.avg, breaks = length(col), include.lowest = FALSE)]
# 
#   return(colors)
# }
# 
# # example of estimation
# fLLMgn<-TVAMLLMgnReg(Y, t, x, T, X, h0=h0, h=h)
# 
# # example plot for marginal f0
# par(mfrow=c(1,1))
# plot(t,cos(2*pi*t)*0.5,type='l',lwd=2,lty=4)
# points(t,fLLMgn$f0[,1],type='l',lwd=2,col=2)
# abline(h=0,col=8)
# 
# # example plots for marginal fj
# f1<-f2<-matrix(nrow=M,ncol=N)
# for(m in 1:M){
#   for(k in 1:N){
#     f1[m,k]<-cos(2*pi*t[m])*x0[k]
#   }
#   f2[m,]<-sum(f1[m,])*diff(x0)[1]
# }
# 
# par(mfrow=c(2,2))
# persp(t,x0,f1,theta=-30,phi=30, col = surf.colors(f1), lty=2, border=NA, shade=0.15)
# persp(t,x0,f2,theta=-30,phi=30, col = surf.colors(f2), lty=2, border=NA, shade=0.15)
# persp(t,x0,fLLMgn$f[,,1,1],theta=-30,phi=30, col = surf.colors(fLLMgn$f[,,1,1]), lty=2, border=NA, shade=0.15)
# persp(t,x0,fLLMgn$f[,,1,2],theta=-30,phi=30, col = surf.colors(fLLMgn$f[,,1,2]), lty=2, border=NA, shade=0.15)
# 


