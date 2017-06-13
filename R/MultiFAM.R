#' Iterative Smooth Backfitting Algorithm
#'
#' Smooth backfitting procedure for functional additive models with multiple predictor processes
#'
#' @param Y An \emph{n}-dimensional vector whose elements consist of scalar responses.
#' @param Xi A \emph{d}-dimentional list whose components consist of two lists of \emph{n} vectors containing the obervation time and functional predictor values. See \code{FPCA} for detail.
#' @param h A \emph{d}-dimensional vector of bandwidths for kernel smoothing to estimate each component function. Common bandwidths are applied for component functions corresponding to the same processes, respectively.
#' @param K A \code{function} object representing the kernel to be used in the smooth backfitting (default is 'epan', the the Epanechnikov kernel.).
#' @param supp A \emph{d} by 2 matrix whose row vectors consist of the lower and upper limits of estimation intervals for each component function (default is the \emph{d}-dimensional rectangle of \emph{[-2,2]}).
#' @param FVE A \emph{d}-dimesional vector whose components consist of the fraction of variance explained of of FPCA for each predictor process. (default is 0.95.)
#'
#' @details \code{MultiFAM} fits functional additive models for a scalar response and multiple predictor processes based on the smooth backfitting algorithm proposed by Han et al. (2016) that \deqn{E(Y | \mathbf{X}) = \sum_{j=1}^d \sum_{k=1}^{K_j} g_{jk}(\xi_{jk}),} where \eqn{\xi_{jk}} stand for the k-th FPC score of the j-th predictor process. \code{MultiFAM} only focuses on mutiple predictor processes case. In fact, the case of univariate predictor is the same with functional additive model proposed by Mueller and Yao (2008). Especially in this development, one can designate an estimation support of additive models when the additive modeling is only allowed over restricted intervals or one is interested in the modeling over the support (see Han et al., 2016).
#'
#' @return A list containing the following fields:
#' \item{xi}{An \emph{N} by \emph{d} matrix whose column vectors consist of FPC score grid vectors that each additive component functional is evluated.}
#' \item{SBFit}{An \emph{N} by \emph{d} matrix whose column vectors consist of the smooth backfitting component function estimators at the given estimation points.}
#' \item{mY}{A scalar of centered part of the regression model.}
#' \item{NW}{An \emph{N} by \emph{d} matrix whose column vectors consist of the Nadaraya-Watson marginal regression function estimators for each predictor component at the given estimation points.}
#' \item{mgnDens}{An \emph{N} by \emph{d} matrix whose column vectors consist of the marginal kernel density estimators for each predictor component at the given estimation points.}
#' \item{jntDens}{An \emph{N} by \emph{N} by \emph{d} by \emph{d} array representing the 2-dimensional joint kernel density estimators for all pairs of predictor components at the given estimation grid. For example, \code{[,,j,k]} of the object provides the 2-dimensional joint kernel density estimator of the \code{(j,k)}-component of predictor components at the corresponding \emph{N} by \emph{N} matrix of estimation grid.}
#' \item{itemNum}{The iteration number that the smooth backfitting algorithm has stopped.}
#' \item{itemErr}{The iteration error of the smooth backfitting algorithm that represents the maximum L2 distance among component functions in the last successive updates.}
#' @examples
#' set.seed(1000)
#' 
#' library(MASS)
#' 
#' f11 <- function(t) t
#' f12 <- function(t) 4*sin(3/2*pi*t/4)
#' f21 <- function(t) 3*atan(2*pi*t/4)
#' f22 <- function(t) 4*cos(2*pi*t/4)
#' 
#' n<-250
#' 
#' sig <- matrix(c(1.5, 0.0, 0.7, -.2,
#'                 0.0, 1.2, -.3, 0.6,
#'                 0.7, -.3, 2.0, 0.0,
#'                 -.2, 0.6, 0.0, 1.5),
#'               nrow=4,ncol=4)
#' 
#' scoreX <- mvrnorm(n,mu=rep(0,4),Sigma=sig)
#' Y <- f11(scoreX[,1]) + f12(scoreX[,2]) + f21(scoreX[,3]) + f22(scoreX[,4]) + rnorm(n,0,0.5)
#' 
#' phi11 <- function(t) sqrt(2)*sin(2*pi*t)
#' phi12 <- function(t) sqrt(2)*sin(4*pi*t)
#' phi21 <- function(t) sqrt(2)*cos(2*pi*t)
#' phi22 <- function(t) sqrt(2)*cos(4*pi*t)
#' 
#' grid <- seq(0,1,length.out=51)
#' Lt <- Lx1 <- Lx2 <- list()
#' for (i in 1:n) {
#'   Lt[[i]] <- grid
#'   Lx1[[i]] <- scoreX[i,1]*phi11(grid) + scoreX[i,2]*phi12(grid)
#'   Lx2[[i]] <- scoreX[i,3]*phi21(grid) + scoreX[i,4]*phi22(grid)
#' }
#' 
#' Xi1 <- list(Ly=Lx1, Lt=Lt)
#' Xi2 <- list(Ly=Lx2, Lt=Lt)
#' 
#' Xi <- list(Xi1, Xi2)
#' 
#' h <- c(0.5, 0.5)
#' mFAMfit <- MultiFAM(Y,Xi,h=h,FVE=NULL)
#' 
#' par(mfrow=c(2,2))
#' plot(mFAMfit$xi[,1],f11(mFAMfit$xi[,1]),type='l',lty=4,
#'      xlab=expression(xi[11]),ylab='f11',ylim=c(-5,5))
#' points(mFAMfit$xi[,1],-mFAMfit$SBFit[,1],type='l',col=2,lwd=2)
#' abline(h=0,col=8)
#' 
#' plot(mFAMfit$xi[,2],f12(mFAMfit$xi[,2]),type='l',lty=4,
#'      xlab=expression(xi[12]),ylab='f12',ylim=c(-5,5))
#' points(mFAMfit$xi[,2],mFAMfit$SBFit[,2],type='l',col=2,lwd=2)
#' abline(h=0,col=8)
#' 
#' plot(mFAMfit$xi[,3],f21(mFAMfit$xi[,3]),type='l',lty=4,
#'      xlab=expression(xi[21]),ylab='f21',ylim=c(-5,5))
#' points(mFAMfit$xi[,3],mFAMfit$SBFit[,3],type='l',col=2,lwd=2)
#' abline(h=0,col=8)
#' 
#' plot(mFAMfit$xi[,4],f22(mFAMfit$xi[,4]),type='l',lty=4,
#'      xlab=expression(xi[22]),ylab='f22',ylim=c(-5,5))
#' points(mFAMfit$xi[,4],mFAMfit$SBFit[,4],type='l',col=2,lwd=2)
#' abline(h=0,col=8)
#' 
#' @references
#' \cite{Mammen, E., Linton, O. and Nielsen, J. (1999), "The existence and asymptotic properties of a backfitting projection algorithm under weak conditions", Annals of Statistics, Vol.27, No.5, p.1443-1490.}
#'
#' \cite{Mammen, E. and Park, B. U. (2006), "A simple smooth backfitting method for additive models", Annals of Statistics, Vol.34, No.5, p.2252-2271.}
#'
#' \cite{Mueller, H.-G. and Yao, F. (2008), "Functional additive models", Journal of the Americal Statistical Association, Vol.103, No.484, p.1534-1544.}
#'
#' \cite{Han, K., Mueller, H.-G. and Park, B. U. (2016), "Smooth backfitting for additive modeling with small errors-in-variables, with an application to additive functional regression for multiple predictor functions", Bernoulli (accepted).}
#' @export

MultiFAM <- function(Y,Xi,h=NULL,K='epan',FVE=NULL,supp=NULL){

  if (length(Xi)<2) {
    return(message('Predictor process must be multi-dimensional.'))
  }
  
  N <- 51
  n <- length(Y)
  d0 <- length(Xi)
  
  if (is.null(FVE)==TRUE) {
    FVE <- rep(0.95,d0)
  } else {
    if (length(FVE)<2) {
      FVE <- rep(FVE,d0)
    }
  }
  if (is.null(supp)==TRUE) {
    supp <- matrix(rep(c(-2,2),d0),ncol=2,byrow=TRUE)
  }
  
  x <- c()
  
  dj <- c()
  X <- c()
  xi <- c()
  supp0 <- matrix(ncol=2)
  for (j in 1:d0) {
    tmpLy <- Xi[[j]]$Ly
    tmpLt <- Xi[[j]]$Lt
    
    tmpFPCA <- FPCA(tmpLy, tmpLt, optns = list(FVEthreshold=FVE[j]))
    
    Xij <- t(t(tmpFPCA$xiEst)/sqrt(tmpFPCA$lambda))
    dj <- length(tmpFPCA$lambda)
    
    X <- cbind(X,Xij)
    
    xj0 <- seq(supp[j,1],supp[j,2],length.out=N)
    x <- cbind(x, matrix(rep(xj0,dj),ncol=dj))
    
    xi0 <- matrix(rep(xj0,dj),nrow=N,ncol=dj)
    xij <- t(t(xi0)*sqrt(tmpFPCA$lambda))
    xi <- cbind(xi,xij)
    
    supp0 <- t(cbind(t(supp0),matrix(rep(supp[j,],dj),ncol=2)))
  }
  
  supp <- supp0[-1,]
  
  d <- ncol(X)
  
  if (K!='epan') {
    message('Epanechnikov kernel is only supported currently. It uses Epanechnikov kernel automatically')
    K<-'epan'
  }
  if (is.null(h)==TRUE) {
    h <- rep(0.25*n^(-1/5),d)*(supp[,2]-supp[,1])
  } else {
    h <- rep(h,dj)
  }
  if (length(h)<2) {
    return(message('Bandwidth must be multi-dimensional.'))
  }
  
  tmpIndex <- rep(1,n)
  for (l in 1:d) {
    tmpIndex <- tmpIndex*dunif(X[,l],supp[l,1],supp[l,2])*(supp[l,2]-supp[l,1])
  }
  tmpIndex <- which(tmpIndex==1)
  
  yMean <- sum(Y[tmpIndex])/length(Y)/P0(X,supp)   
  
  MgnJntDens <- MgnJntDensity(x,X,h,K,supp)
  
  fNW <- NWMgnReg(Y, x, X, h, K, supp)

  # for (j in 1:d) {
  #   plot(x[,j],fNW[,j],type='l')
  # }
  
  f <- matrix(0,nrow=N,ncol=d)
  
  # backfitting
  eps <- epsTmp <- 100
  iter <- 1
  
  critEps <- 5e-5
  critEpsDiff <- 5e-4
  critIter <- 100
  
  while (eps>critEps) {
    #print(eps)
    epsTmp <- eps
    f0 <- f
    
    for (j in 1:d) {
      f[,j] <- SBFCompUpdate(f,j,fNW,Y,X,x,h,K,supp,MgnJntDens)[,j]
      
      if (sum(is.nan(f[,j])==TRUE)>0) {
        f[which(is.nan(f[,j])==TRUE),j] <- 0
      }
      
      
      if (sum(f[,j]*f0[,j])<0) {
        f[,j] <- -f[,j]
      }
    }
    
    # for (j in 1:d) {
    #   plot(x[,j],f[,j],type='l')
    # }
    
    eps <- max(sqrt(apply(abs(f-f0)^2,2,'mean')))
    
    if (abs(epsTmp-eps)<critEpsDiff) {
      return(list(
        xi=xi,
        SBFit=f, 
        mY=yMean,
        NW=fNW, 
        mgnDens=MgnJntDens$pMatMgn, 
        jntDens=MgnJntDens$pArrJnt, 
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
        xi=xi,
        SBFit=f, 
        mY=yMean,
        NW=fNW, 
        mgnDens=MgnJntDens$pMatMgn, 
        jntDens=MgnJntDens$pArrJnt, 
        iterNum=iter,
        iterErr=eps,
        iterErrDiff=abs(epsTmp-eps),
        critNum=critIter,
        critErr=critEps,
        critErrDiff=critEpsDiff
      )
      )
    }
    
    iter <- iter+1
  }
  
  return(list(
    xi=xi,
    SBFit=f, 
    mY=yMean,
    NW=fNW, 
    mgnDens=MgnJntDens$pMatMgn, 
    jntDens=MgnJntDens$pArrJnt, 
    iterNum=iter,
    iterErr=eps,
    iterErrDiff=abs(epsTmp-eps),
    critNum=critIter,
    critErr=critEps,
    critErrDiff=critEpsDiff
  )
  )
  
}


