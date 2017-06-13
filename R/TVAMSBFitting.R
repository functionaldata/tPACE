#' Iterative Smooth Backfitting Algorithm
#'
#' Smooth backfitting procedure for nonparametric time-varying additive models
#'
#' @param Y An \emph{n}-dimensional list whose elements consist of vectors of longitudinal responses.
#' @param t An \emph{M}-dimensional vector of estimation time points for marginal mean curve of \emph{Y(t)}.
#' @param x An \emph{N} by \emph{d} matrix whose column vectors consist of \emph{N} vectors of estimation points for each component surface.
#' @param T An \emph{n}-dimensional list whose elements consist of vectors of time points corresponding to longitudinal responses \emph{Y}.
#' @param X An \emph{n} by \emph{d} matrix whose row vectors consist of multivariate predictors.
#' @param h0 A scalar value of bandwidth for kernel smoothing to estimate marginal mean curve of \emph{Y(t)}.
#' @param h A \emph{d}-dimensional vector of bandwidths for kernel smoothing to estimate each component surface.
#' @param K A \code{function} object representing the kernel to be used in the smooth backfitting (default is 'epan', the the Epanechnikov kernel.).
#' @param supp0 A 2-dimensional vector whose elements consist of the lower and upper limit of estimation interval for marginal mean curve (default is the \emph{1}-dimensional unit interval, \emph{[0,1]}).
#' @param supp A \emph{d} by 2 matrix whose row vectors consist of the lower and upper limits of estimation intervals for each component function (default is the \emph{d}-dimensional unit rectangle of \emph{[0,1]^d}).
#'
#' @details \code{TVAMSBFitting} fits component surface of time-varying additive models for a longitudinal response and a multivariate predictor based on the smooth backfitting algorithm, a simplified model of Zhang et al. (2013). Here, the time-varying additive model stands for \deqn{E(Y(t) | \mathbf{X}) = g_1(t,X_1)+\cdots+ g_d(t,X_d).} \code{TVAMSBFitting} only focuses on the local linear smooth backfitting estimator with multivariate predictor case for the purpose of enjoying smooth surface estimation. Support of the multivariate predictor is assumed to be a product of closed intervals. Especially in this development, one can designate an estimation support of additive models when the additive modeling is only allowed over restricted intervals or one is interested in the modeling over the support (see Han et al., 2016). If one puts \code{X} on the argument of estimation points \code{x}, \code{TVAMSBFitting} returns estimated values of conditional mean curve for responses given observed predictors.
#'
#' @return A list containing the following fields:
#' \item{SBFit}{An \emph{N} by \emph{N} by \emph{d} array whose elements consist of the locally linear smooth backfitting component surface estimators at the given estimation points.}
#' \item{f0}{An \emph{M} vector whose elements consist of marginal mean curve of \emph{Y(t)}.}
#' \item{fj}{An \emph{N} by \emph{N} by \emph{d} array whose elements consist of the locally linear kernel smoothing estimators of marginal component surface.}
#' \item{itemNum}{The iteration number that the smooth backfitting algorithm has stopped.}
#' \item{itemErr}{The iteration error of the smooth backfitting algorithm that represents the maximum L2 distance among component functions in the last successive updates.}
#' @examples
#' set.seed(100)
#' 
#' n <- 200
#' N <- 51
#' M <- 81
#' d <- 2
#' 
#' nSet <- sample(50:100,n,replace=TRUE)
#' T <- list()
#' for (i in 1:n){
#'   T[[i]] <- runif(nSet[i],0,1)
#' }
#' 
#' X <- pnorm(matrix(rnorm(n*d),nrow=n,ncol=d)%*%matrix(c(1,0.5,0.5,1),nrow=2,ncol=2))
#' 
#' g1 <- function(x1) 2*(x1-0.5)
#' g2 <- function(x2) sin(2*pi*x2)
#' 
#' g <- function(u,x) sin(2*pi*u)*g1(x[1]) + cos(2*pi*u)*g2(x[2])
#' 
#' Y <- list()
#' for (i in 1:n){
#'   tmpY <- c()
#'   for (j in 1:nSet[i]){
#'     tmpY[j] <- g(T[[i]][j],X[i,])
#'   }
#'   Y[[i]] <- tmpY+rnorm(nSet[i],0,0.2)
#' }
#' 
#' t <- seq(0,1,length.out=M)
#' x <- matrix(rep(seq(0,1,length.out=N),d),nrow=N,ncol=d)
#' 
#' h0 <- 0.15
#' h <- c(0.15,0.15)
#' 
#' sbfLLSurf <- TVAMSBFitting(Y,t,x,T,X,h0=h0,h=h)
#' 
#' g1SbfLL <- sbfLLSurf$SBFit[,,1]
#' g2SbfLL <- sbfLLSurf$SBFit[,,2]
#' 
#' g1Eval <- matrix(nrow=M,ncol=N)
#' g2Eval <- matrix(nrow=M,ncol=N)
#' for (i in 1:M) {
#'   for (j in 1:N) {
#'     g1Eval[i,j] <- sin(2*pi*t[i])*g1(x[j,1])
#'     g2Eval[i,j] <- cos(2*pi*t[i])*g2(x[j,2])
#'   }
#' }
#' 
#' par(mfrow=c(2,2))
#' persp(t,x[,1],g1Eval,theta=-30,phi=30,lty=2, border=NA, shade=0.5,xlab='t',ylab='x1',zlab='g1Eval')
#' persp(t,x[,2],g2Eval,theta=20,phi=30,lty=2, border=NA, shade=0.5,xlab='t',ylab='x2',zlab='g2Eval')
#' 
#' persp(t,x[,1],g1SbfLL,theta=-30,phi=40, lty=2, border=NA, shade=0.5,xlab='t',ylab='x1',zlab='g1Est')
#' persp(t,x[,2],g2SbfLL,theta=20,phi=40, lty=2, border=NA, shade=0.5,xlab='t',ylab='x2',zlab='g2Est')
#' @references
#' \cite{Mammen, E., Linton, O. and Nielsen, J. (1999), "The existence and asymptotic properties of a backfitting projection algorithm under weak conditions", Annals of Statistics, Vol.27, No.5, p.1443-1490.}
#'
#' \cite{Mammen, E. and Park, B. U. (2006), "A simple smooth backfitting method for additive models", Annals of Statistics, Vol.34, No.5, p.2252-2271.}
#'
#' \cite{Yu, K., Park, B. U. and Mammen, E. (2008), "Smooth backfitting in generalized additive models", Annals of Statistics, Vol.36, No.1, p.228-260.}
#'
#' \cite{Lee, Y. K., Mammen, E. and Park., B. U. (2010), "backfitting and smooth backfitting for additive quantile models", Annals of Statistics, Vol.38, No.5, p.2857-2883.}
#'
#' \cite{Lee, Y. K., Mammen, E. and Park., B. U. (2012), "Flexible generalized varying coefficient regression models", Annals of Statistics, Vol.40, No.3, p.1906-1933.}
#'
#' \cite{Han, K., Mueller, H.-G. and Park, B. U. (2016), "Smooth backfitting for additive modeling with small errors-in-variables, with an application to additive functional regression for multiple predictor functions", Bernoulli (accepted).}
#'
#' \cite{Zhang, X., Park, B. U. and Wang, J.-L. (2013), "Time-varying additive models for longitudinal data", Journal of the American Statistical Association, Vol.109, No.503, p.983-998.}
#' @export

TVAMSBFitting <- function(Y,t,x,T,X,h0=NULL,h=NULL,K='epan',supp0=NULL,supp=NULL){
  
  if (is.null(ncol(x))==TRUE) {
    return(message('Evaluation grid must be multi-dimensional.'))
  }
  if (is.null(ncol(X))==TRUE) {
    return(message('Observation grid must be multi-dimensional.'))
  }
  if (length(h)<2) {
    return(message('Bandwidth must be multi-dimensional.'))
  }
  
  N <- nrow(x)
  d <- ncol(x)
  n <- length(Y)
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
  
  MgnJntDens <- TVAMMgnJntDensity(t, x, T, X, h0=h0, h=h)
  
  fLL <- TVAMLLMgnReg(Y, t, x, T, X, h0=h0, h=h)
  f0LL <- fLL$f0  
  fjLL <- fLL$f

  f <- array(0,c(M,N,3,d))
  
  # backfitting
  eps <- epsTmp <- 100
  iter <- 1
  
  critEps <- 5e-4
  critEpsDiff <- 5e-3
  critIter <- 100
  
  while (eps>critEps) {
    #print(eps)
    
    epsTmp <- eps
    f0 <- f
    
    for (j in 1:d) {
      f[,,,j] <- TVAMSBFCompUpdate(f,j,fLL,Y,T,X,t,x,h0,h,K,supp0,supp,MgnJntDens)[,,,j]
      
      if (sum(is.nan(f[,,,j])==TRUE)>0) {
        f[which(is.nan(f[,,,j])==TRUE,arr.ind=TRUE),j] <- 0
      }
      
      if (sum(f[,,1,j]*f0[,,1,j])<0) {
        f[,,,j] <- -f[,,,j]
      }
    }
    
    # eps0 <- max(apply(abs(f[,,1,]-f0[,,1,]),3,'mean'))
    # eps1 <- max(apply(abs(f[,,2,]-f0[,,2,]),3,'mean'))/h0
    # eps2 <- max(apply(abs(f[,,3,]-f0[,,3,]),3,'mean'))/min(h)
    # 
    # #print(c(eps0,eps1,eps2))
    # 
    # tmp1 <- max(c(eps0,eps1,eps2))
    # tmp2 <- min(c(eps0,eps1,eps2))
    # 
    # print(c(tmp1,tmp2))
    # 
    # eps <- mean(c(tmp1,tmp2))
    
    #eps <- max(abs(f[,,1,]-f0[,,1,]))
    eps <- max(sqrt(apply(abs(f[,,1,]-f0[,,1,])^2,2,'mean')))
        
    if (abs(epsTmp-eps)<critEpsDiff) {
      return(list(
          SBFit=f[,,1,], 
          f0=f0LL[,1],
          fj=fjLL[,,1,],
          mgnDens=MgnJntDens$qArrMgn, 
          jntDens=MgnJntDens$qArrJnt, 
          iterNum=iter,
          iterErr=eps,
          iterErrDiff=epsTmp-eps,
          critNum=critIter,
          critErr=critEps,
          critErrDiff=critEpsDiff
        )
      )
    }

    if (iter>critIter) {
      message('The algorithm may not converge (SBF iteration > stopping criterion). Try another choice of bandwidths.')
      return(list(
          SBFit=f[,,1,], 
          f0=f0LL[,1],
          fj=fjLL[,,1,], 
          mgnDens=MgnJntDens$qArrMgn, 
          jntDens=MgnJntDens$qArrJnt, 
          iterNum=iter,
          iterErr=eps,
          iterErrDiff=epsTmp-eps,
          critNum=critIter,
          critErr=critEps,
          critErrDiff=critEpsDiff
        )
      )
    }
    
    iter <- iter+1
  }
  
  return(list(
      SBFit=f[,,1,], 
      f0=f0LL[,1],
      fj=fjLL[,,1,],
      mgnDens=MgnJntDens$qArrMgn, 
      jntDens=MgnJntDens$qArrJnt, 
      iterNum=iter,
      iterErr=eps,
      iterErrDiff=epsTmp-eps,
      critNum=critIter,
      critErr=critEps,
      critErrDiff=critEpsDiff
    )
  )
  
}





# 
# # test example
# M<-81
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
# time1<-Sys.time()
# fit<-TVAMSBFitting(Y,t,x,T,X,h0=h0,h=h)
# time2<-Sys.time()
# 
# time2-time1
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
# f1<-sin(2*pi*t)%*%t(g1(x0))
# f2<-cos(2*pi*t)%*%t(g2(x0))
# f1Fit<-fit$SBFit[,,1,1]
# f2Fit<-fit$SBFit[,,1,2]
# 
# # example of centering
# par(mfrow=c(2,2))
# persp(t,x0,f1,theta=30,phi=40, col = surf.colors(f1), lty=2, border=NA, shade=0.1)
# persp(t,x0,f2,theta=30,phi=40, col = surf.colors(f2), lty=2, border=NA, shade=0.1)
# persp(t,x0,f1Fit,theta=30,phi=40, col = surf.colors(f1Fit), lty=2, border=NA, shade=0.1)
# persp(t,x0,f2Fit,theta=30,phi=40, col = surf.colors(f2Fit), lty=2, border=NA, shade=0.1)