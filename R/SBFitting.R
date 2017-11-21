#' Iterative Smooth Backfitting Algorithm
#'
#' Smooth backfitting procedure for nonparametric additive models
#'
#' @param Y An \emph{n}-dimensional vector whose elements consist of scalar responses.
#' @param x An \emph{N} by \emph{d} matrix whose column vectors consist of \emph{N} vectors of estimation points for each component function.
#' @param X An \emph{n} by \emph{d} matrix whose row vectors consist of multivariate predictors.
#' @param h A \emph{d} vector of bandwidths for kernel smoothing to estimate each component function.
#' @param K A \code{function} object representing the kernel to be used in the smooth backfitting (default is 'epan', the the Epanechnikov kernel.).
#' @param supp A \emph{d} by 2 matrix whose row vectors consist of the lower and upper limits of estimation intervals for each component function (default is the \emph{d}-dimensional unit rectangle of \emph{[0,1]}).
#'
#' @details \code{SBFitting} fits component functions of additive models for a scalar response and a multivariate predictor based on the smooth backfitting algorithm proposed by Mammen et al. (1999) and intensively studied by Mammen and Park (2006), Yu et al. (2008), Lee et al. (2010, 2012) and so on. \code{SBFitting} only focuses on the local constant smooth backfitting estimator with multivariate predictor case. In fact, the case of univariate predictor is the same with the local constant kernel regression estimator (Nadaraya-Watson estimator) and the local polynomial version can be extended similarly, so that those are omitted in the development. Support of the multivariate predictor is assumed to be a product of closed intervals. Especially in this development, one can designate an estimation support of additive models when the additive modeling is only allowed over restricted intervals or one is interested in the modeling over the support (see Han et al., 2016). If one puts \code{X} on the argument of estimation points \code{x}, \code{SBFitting} returns estimated values of conditional mean responses given observed predictors.
#'
#' @return A list containing the following fields:
#' \item{SBFit}{An \emph{N} by \emph{d} matrix whose column vectors consist of the smooth backfitting component function estimators at the given estimation points.}
#' \item{mY}{A scalar of centered part of the regression model.}
#' \item{NW}{An \emph{N} by \emph{d} matrix whose column vectors consist of the Nadaraya-Watson marginal regression function estimators for each predictor component at the given estimation points.}
#' \item{mgnDens}{An \emph{N} by \emph{d} matrix whose column vectors consist of the marginal kernel density estimators for each predictor component at the given estimation points.}
#' \item{jntDens}{An \emph{N} by \emph{N} by \emph{d} by \emph{d} array representing the 2-dimensional joint kernel density estimators for all pairs of predictor components at the given estimation grid. For example, \code{[,,j,k]} of the object provides the 2-dimensional joint kernel density estimator of the \code{(j,k)}-component of predictor components at the corresponding \emph{N} by \emph{N} matrix of estimation grid.}
#' \item{itemNum}{The iteration number that the smooth backfitting algorithm has stopped.}
#' \item{itemErr}{The iteration error of the smooth backfitting algorithm that represents the maximum L2 distance among component functions in the last successive updates.}
#' @examples
#' set.seed(100)
#' 
#' n <- 100
#' d <- 2
#' X <- pnorm(matrix(rnorm(n*d),nrow=n,ncol=d)%*%matrix(c(1,0.6,0.6,1),nrow=2,ncol=2))
#' 
#' f1 <- function(t) 2*(t-0.5)
#' f2 <- function(t) sin(2*pi*t)
#' 
#' Y <- f1(X[,1])+f2(X[,2])+rnorm(n,0,0.1)
#' 
#' # component function estimation
#' N <- 101
#' x <- matrix(rep(seq(0,1,length.out=N),d),nrow=N,ncol=d)
#' h <- c(0.12,0.08)
#'   
#' sbfEst <- SBFitting(Y,x,X,h)
#' fFit <- sbfEst$SBFit
#' 
#' par(mfrow=c(1,2))
#' plot(x[,1],f1(x[,1]),type='l',lwd=2,col=2,lty=4,xlab='X1',ylab='Y')
#' points(x[,1],fFit[,1],type='l',lwd=2,col=1)
#' points(X[,1],Y,cex=0.3,col=8)
#' legend('topleft',legend=c('SBF','true'),col=c(1,2),lwd=2,lty=c(1,4),horiz=FALSE,bty='n')
#' abline(h=0,col=8)
#' 
#' plot(x[,2],f2(x[,2]),type='l',lwd=2,col=2,lty=4,xlab='X2',ylab='Y')
#' points(x[,2],fFit[,2],type='l',lwd=2,col=1)
#' points(X[,2],Y,cex=0.3,col=8)
#' legend('topright',legend=c('SBF','true'),col=c(1,2),lwd=2,lty=c(1,4),horiz=FALSE,bty='n')
#' abline(h=0,col=8)
#' 
#' # prediction
#' x <- X
#' h <- c(0.12,0.08)
#'   
#' sbfPred <- SBFitting(Y,X,X,h)
#' fPred <- sbfPred$mY+apply(sbfPred$SBFit,1,'sum')
#' 
#' par(mfrow=c(1,1))
#' plot(fPred,Y,cex=0.5,xlab='SBFitted values',ylab='Y')
#' abline(coef=c(0,1),col=8)
#' @references
#' \cite{Mammen, E., Linton, O. and Nielsen, J. (1999), "The existence and asymptotic properties of a backfitting projection algorithm under weak conditions", Annals of Statistics, Vol.27, No.5, p.1443-1490.}
#'
#' \cite{Mammen, E. and Park, B. U. (2006), "A simple smooth backfitting method for additive models", Annals of Statistics, Vol.34, No.5, p.2252-2271.}
#'
#' \cite{Yu, K., Park, B. U. and Mammen, E. (2008), "Smooth backfitting in generalized additive models", Annals of Statistics, Vol.36, No.1, p.228-260.}
#'
#' \cite{Lee, Y. K., Mammen, E. and Park., B. U. (2010), "backfitting and smooth backfitting for additive quantile models", Vol.38, No.5, p.2857-2883.}
#'
#' \cite{Lee, Y. K., Mammen, E. and Park., B. U. (2012), "Flexible generalized varying coefficient regression models", Annals of Statistics, Vol.40, No.3, p.1906-1933.}
#'
#' \cite{Han, K., Mueller, H.-G. and Park, B. U. (2016), "Smooth backfitting for additive modeling with small errors-in-variables, with an application to additive functional regression for multiple predictor functions", Bernoulli (accepted).}
#' @export

SBFitting <- function(Y,x,X,h=NULL,K='epan',supp=NULL){
  
  if (is.null(ncol(x))==TRUE) {
    return(message('Evaluation grid must be multi-dimensional.'))
  }
  if (is.null(ncol(X))==TRUE) {
    return(message('Observation grid must be multi-dimensional.'))
  }
  
  N <- nrow(x)
  n <- nrow(X)
  d <- ncol(X)
  
  if (K!='epan') {
    message('Epanechnikov kernel is only supported currently. It uses Epanechnikov kernel automatically')
    K<-'epan'
  }
  if (is.null(supp)==TRUE) {
    supp <- matrix(rep(c(0,1),d),ncol=2,byrow=TRUE)
  }
  if (is.null(h)==TRUE) {
    h <- rep(0.25*n^(-1/5),d)*(supp[,2]-supp[,1])
  }
  if (length(h)<2) {
    return(message('Bandwidth must be multi-dimensional.'))
  }
  
  tmpIndex <- rep(1,n)
  for (l in 1:d) {
    tmpIndex <- tmpIndex*dunif(X[,l],supp[l,1],supp[l,2])*(supp[l,2]-supp[l,1])
  }
  tmpIndex <- which(tmpIndex==1)
  
  yMean <- sum(Y[tmpIndex])/length(Y)/P0(X,supp=supp)   
  
  MgnJntDens <- MgnJntDensity(x,X,h,K,supp)
  
  fNW <- NWMgnReg(Y, x, X, h, K, supp)
  f <- matrix(0,nrow=N,ncol=d)
  
  # backfitting
  eps <- epsTmp <- 100
  iter <- 1
  
  critEps <- 5e-5
  critEpsDiff <- 5e-4
  critIter <- 50
  
  while (eps>critEps) {
    
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
    
    #eps <- max(abs(f-f0))
    eps <- max(sqrt(apply(abs(f-f0)^2,2,'mean')))
    
    if (abs(epsTmp-eps)<critEpsDiff) {
      return(list(
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
