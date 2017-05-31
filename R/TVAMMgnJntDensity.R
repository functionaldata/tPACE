#####
##### marginal and 2-dim. joint kernel densities estimators
#####

##### input variables:
#####   j: index of kernel estimation for marginal density (scalar)
#####   kj: index of kernel estimation for 2-dim. joint density (2-dim. vector)
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
#####   q0: margianl densities of time at estimation points near observation points (M-dim. vector)
#####   qArrMgn: Locally linear margianl densities of time and j-th covariate at estimation points near observation points (M*N*3*3*d array)
#####   qArrJnt: Locally linear joint densities of time, k-th and j-th covariate at estimation grid near observation grid (M*N*N*3*3*d*d array)


Q0 <- function(t, T, h0=NULL, K='epan', supp0=NULL){
  
  n <-length(T)
  M <- length(t)
  
  if (K!='epan') {
    message('Epanechnikov kernel is only supported currently. It uses Epanechnikov kernel automatically')
    K<-'epan'
  }
  if (is.null(supp0)==TRUE) {
    supp0 <- c(0,1)
  }
  if (is.null(h0)==TRUE) {
    h0 <- 0.25*n^(-1/5)*(supp0[2]-supp0[1])
  }
  
  vecT <- c()
  for (i in 1:n) {
    vecT<-c(vecT,T[[i]])
  }
  
  normKernelT <- NormKernel(t,vecT,h0,K,supp0)
  qHat0 <-apply(normKernelT,1,'mean')
  
  return(qHat0)
}


### pseudo marginal density for time and one conponent
Qj <- function(j, t, x, T, X, h0=NULL, h=NULL, K='epan', supp0=NULL, supp=NULL){
  
  N <- nrow(x)
  d <- ncol(x)
  n <- nrow(X)
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
  
  vecT <- c()
  tmpRep <- c()
  for (i in 1:n) {
    vecT <- c(vecT,T[[i]])
    tmpRep[i] <- length(T[[i]])
  }
  #vecXj <- rep(X[,j],rep(M,n))
  vecXj <- rep(X[,j],tmpRep)
  
  v0 <- (matrix(rep(vecT,M), nrow=M, ncol=length(vecT), byrow=TRUE)  - matrix(rep(t,length(vecT)), nrow=M, ncol=length(vecT)))/h0
  vj <- (matrix(rep(vecXj,N), nrow=N, ncol=length(vecXj), byrow=TRUE)  - matrix(rep(x[,j],length(vecXj)), nrow=N, ncol=length(vecXj)))/h[j]
  
  normKernelT <- NormKernel(t,vecT,h0,K,supp0)
  normKernelXj <- NormKernel(x[,j],vecXj,h[j],K,supp[j,])
  
  qHatj11 <- normKernelT%*%t(normKernelXj)/length(vecT)
  qHatj12 <- (normKernelT*v0)%*%t(normKernelXj)/length(vecT)
  qHatj13 <- normKernelT%*%t(normKernelXj*vj)/length(vecT)
  qHatj22 <- (normKernelT*v0^2)%*%t(normKernelXj)/length(vecT)
  qHatj23 <- (normKernelT*v0)%*%t(normKernelXj*vj)/length(vecT)
  qHatj33 <- normKernelT%*%t(normKernelXj*vj^2)/length(vecT)
  
  qHatj <- array(cbind(qHatj11,qHatj12,qHatj13,
                       qHatj12,qHatj22,qHatj23,
                       qHatj13,qHatj23,qHatj33),c(M,N,3,3))

  return(qHatj)   
  
}


# pseudo joint density for time and two conponents
Qkj <- function(kj, t, x, T, X, h0=NULL, h=NULL, K='epan', supp0=NULL, supp=NULL){
  
  N <- nrow(x)
  d <- ncol(x)
  n <- nrow(X)
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
  
  vecT <- c()
  tmpRep <- c()
  for (i in 1:n) {
    vecT <- c(vecT,T[[i]])
    tmpRep[i] <- length(T[[i]])
  }

  k <- kj[1]
  #vecXk <- rep(X[,k],rep(M,n))
  vecXk <- rep(X[,k],tmpRep)
  
  j <- kj[2]
  #vecXj <- rep(X[,j],rep(M,n))
  vecXj <- rep(X[,j],tmpRep)
  
  v0 <- (matrix(rep(vecT,M), nrow=M, ncol=length(vecT), byrow=TRUE)  - matrix(rep(t,length(vecT)), nrow=M, ncol=length(vecT)))/h0
  vk <- (matrix(rep(vecXk,N), nrow=N, ncol=length(vecXk), byrow=TRUE)  - matrix(rep(x[,k],length(vecXk)), nrow=N, ncol=length(vecXk)))/h[k]
  vj <- (matrix(rep(vecXj,N), nrow=N, ncol=length(vecXj), byrow=TRUE)  - matrix(rep(x[,j],length(vecXj)), nrow=N, ncol=length(vecXj)))/h[j]
  
  normKernelT <- NormKernel(t,vecT,h0,K,supp0)
  normKernelXk <- NormKernel(x[,k],vecXk,h[k],K,supp[k,])
  normKernelXj <- NormKernel(x[,j],vecXj,h[j],K,supp[j,])
  
  ind1 <- max(which(t<2*h0))
  ind2 <- min(which(t>1-2*h0))

  qHatkj <- array(0,dim=c(M,N,N,3,3))
  for (m in 1:M) {
    #print(m)
    
    qHatkj11 <- normKernelXk%*%Matrix::sparseMatrix(1:length(vecT),1:length(vecT),x=c(normKernelT[m,]))%*%t(normKernelXj)/length(vecT)
    qHatkj12 <- normKernelXk%*%Matrix::sparseMatrix(1:length(vecT),1:length(vecT),x=c(normKernelT[m,]*v0[m,]))%*%t(normKernelXj)/length(vecT)
    qHatkj13 <- normKernelXk%*%Matrix::sparseMatrix(1:length(vecT),1:length(vecT),x=c(normKernelT[m,]))%*%t(normKernelXj*vj)/length(vecT)
    
    qHatkj21 <- normKernelXk%*%Matrix::sparseMatrix(1:length(vecT),1:length(vecT),x=c(normKernelT[m,]*v0[m,]))%*%t(normKernelXj)/length(vecT)
    qHatkj22 <- normKernelXk%*%Matrix::sparseMatrix(1:length(vecT),1:length(vecT),x=c(normKernelT[m,]*v0[m,]^2))%*%t(normKernelXj)/length(vecT)
    qHatkj23 <- normKernelXk%*%Matrix::sparseMatrix(1:length(vecT),1:length(vecT),x=c(normKernelT[m,]*v0[m,]))%*%t(normKernelXj*vj)/length(vecT)
    
    qHatkj31 <- (normKernelXk*vk)%*%Matrix::sparseMatrix(1:length(vecT),1:length(vecT),x=c(normKernelT[m,]))%*%t(normKernelXj)/length(vecT)
    qHatkj32 <- (normKernelXk*vk)%*%Matrix::sparseMatrix(1:length(vecT),1:length(vecT),x=c(normKernelT[m,]*v0[m,]))%*%t(normKernelXj)/length(vecT)
    qHatkj33 <- (normKernelXk*vk)%*%Matrix::sparseMatrix(1:length(vecT),1:length(vecT),x=c(normKernelT[m,]))%*%t(normKernelXj*vj)/length(vecT)
    
    qHatkj[m,,,,] <- array(cbind(qHatkj11,qHatkj21,qHatkj31,
                                 qHatkj12,qHatkj22,qHatkj32,
                                 qHatkj13,qHatkj23,qHatkj33),c(N,N,3,3))
    
  }
  
  return(qHatkj)     
  
}


# construction of evaluation matrices for marginal and joint densities estimators
TVAMMgnJntDensity <- function(t, x, T, X, h0=NULL, h=NULL, K='epan', supp0=NULL, supp=NULL){
  
  N <- nrow(x)
  d <- ncol(x)
  n <- nrow(X)
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
  
  q0 <- Q0(t,T,h0,K,supp0)
  qArrMgn <- array(0,dim=c(M,N,3,3,d))
  qArrJnt <- array(0,dim=c(M,N,N,3,3,d,d))
  #cat(paste('Computing all pairs of 1-/2-dim.l marginal/joint density estimators...','\n',sep=''))
  for (j in 1:d) {
    #cat(paste('   ',round(j/d,3),'\n',sep=''))
    #cat('\n')
    qArrMgn[,,,,j] <- Qj(j,t,x,T,X,h0,h,K,supp0,supp)
    
    for (k in j:d) { 
      #print(k)
      if (k==j) {
        next
      } else {
        qArrJnt[,,,,,k,j] <- Qkj(c(k,j),t,x,T,X,h0,h,K,supp0,supp)
        qArrJnt[,,,,,j,k] <- aperm(qArrJnt[,,,,,k,j],c(1,3,2,5,4))
        #qArrJnt[,,,,,j,k] <- Qkj(c(j,k),t,x,T,X,h0,h,K,supp0,supp)
      }
    }
  }
  
  return(list(qArrJnt=qArrJnt, qArrMgn=qArrMgn, q0=q0))
  
}





# # test example
# M<-51
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
# h0<-0.05
# h<-c(0.15,0.15)
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
# # # example of marginal estimation
# # qHat1<-Qj(1, t, x, T, X, h0=h0, h=h)
# # qHat2<-Qj(2, t, x, T, X, h0=h0, h=h)
# #
# # # example plots for marginal qHatj
# # par(mfrow=c(2,2))
# # plot(density(X[,1],bw=h[1],kernel='epanechnikov'))
# # plot(density(X[,2],bw=h[1],kernel='epanechnikov'))
# # persp(t,x0,qHat1[,,1,1],theta=-30,phi=30, col = surf.colors(qHat1[,,1,1]), lty=2, border=NA, shade=0.15)
# # persp(t,x0,qHat2[,,1,1],theta=-30,phi=30, col = surf.colors(qHat2[,,1,1]), lty=2, border=NA, shade=0.15)
# #
# # time1<-Sys.time()
# # # example of joint estimation
# # qHat12<-Qkj(c(1,2), t, x, T, X, h0=h0, h=h)
# # qHat21<-Qkj(c(2,1), t, x, T, X, h0=h0, h=h)
# # time2<-Sys.time()
# #
# # time2-time1
# #
# # # example plots for marginal qHatkj
# # par(mfrow=c(2,3))
# # persp(x0,x0,qHat12[1,,,1,1],theta=-30,phi=30, col = surf.colors(qHat12[7,,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# # persp(x0,x0,qHat12[length(t)/2,,,1,1],theta=-30,phi=30, col = surf.colors(qHat12[length(t)/2,,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# # persp(x0,x0,qHat12[length(t),,,1,1],theta=-30,phi=30, col = surf.colors(qHat12[length(t),,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# # persp(x0,x0,qHat21[1,,,1,1],theta=-30,phi=30, col = surf.colors(qHat21[7,,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# # persp(x0,x0,qHat21[length(t)/2,,,1,1],theta=-30,phi=30, col = surf.colors(qHat21[length(t)/2,,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# # persp(x0,x0,qHat21[length(t),,,1,1],theta=-30,phi=30, col = surf.colors(qHat21[length(t),,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# 
# time1<-Sys.time()
# asdf<-TVAMMgnJntDensity(t, x, T, X, h0=NULL, h=NULL)
# time2<-Sys.time()
# 
# time2-time1
# 
# qHat0<-asdf$q0
# qHat1<-asdf$qArrMgn[,,,,1]
# qHat2<-asdf$qArrMgn[,,,,2]
# 
# par(mfrow=c(1,1))
# plot(t,qHat0,type='l')
# 
# par(mfrow=c(1,2))
# plot(density(X[,1],bw=h[1],kernel='epanechnikov'))
# plot(density(X[,2],bw=h[2],kernel='epanechnikov'))
# persp(t,x0,qHat1[,,1,1],theta=100,phi=30, col = surf.colors(qHat1[,,1,1]), lty=2, border=NA, shade=0.15, ticktype='detailed')
# persp(t,x0,qHat2[,,1,1],theta=100,phi=30, col = surf.colors(qHat2[,,1,1]), lty=2, border=NA, shade=0.15, ticktype='detailed')
# 
# qHat12<-asdf$qArrJnt[,,,,,1,2]
# qHat21<-asdf$qArrJnt[,,,,,2,1]
# 
# par(mfrow=c(2,3))
# persp(x0,x0,qHat12[1,,,1,1],theta=-30,phi=30, col = surf.colors(qHat12[7,,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# persp(x0,x0,qHat12[length(t)/2,,,1,1],theta=-30,phi=30, col = surf.colors(qHat12[length(t)/2,,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# persp(x0,x0,qHat12[length(t),,,1,1],theta=-30,phi=30, col = surf.colors(qHat12[length(t),,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# persp(x0,x0,qHat21[1,,,1,1],theta=-30,phi=30, col = surf.colors(qHat21[7,,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# persp(x0,x0,qHat21[length(t)/2,,,1,1],theta=-30,phi=30, col = surf.colors(qHat21[length(t)/2,,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')
# persp(x0,x0,qHat21[length(t),,,1,1],theta=-30,phi=30, col = surf.colors(qHat21[length(t),,,1,1]), lty=2, border=NA, shade=0.1, ticktype='detailed')


