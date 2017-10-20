#' Functional Additive Models with Multiple Predictor Processes
#'
#' Smooth backfitting procedure for functional additive models with multiple predictor processes
#'
#' @param Y An \emph{n}-dimensional vector whose elements consist of scalar responses.
#' @param Xi A \eqn{(K_1 + ... + K_d)}-dimensional list whose components consist of two lists of \emph{n} vectors containing the obervation time and functional predictor values. See \code{FPCA} for detail.
#' @param K A \code{function} object representing the kernel to be used in the smooth backfitting (default is 'epan', the the Epanechnikov kernel).
#' @param FVE A \emph{d}-dimesional vector whose components consist of the fraction of variance explained of of FPCA for each predictor process (default is 0.95).
#' @param xi An \emph{N} by \eqn{(K_1 + ... + K_d)} matrix whose column vectors consist of \emph{N} vectors of estimation points for each component function.
#' @param h A \eqn{(K_1 + ... + K_d)}-dimensional vector of bandwidths for kernel smoothing to estimate each component function (default is NULL, but it automatically apply shrinkage factor bandwidth selector. See Han et al. (2016)).
#' @param supp A \eqn{(K_1 + ... + K_d)} by 2 matrix whose row vectors consist of the lower and upper limits of estimation intervals for each component function (default is the \emph{d}-dimensional rectangle of \emph{[-2,2]}).
#'
#' @details \code{MultiFAM} fits functional additive models for a scalar response and multiple predictor processes based on the smooth backfitting algorithm proposed by Han et al. (2016) that \deqn{E(Y | \mathbf{X}) = \sum_{j=1}^d \sum_{k=1}^{K_j} g_{jk}(\xi_{jk}),} where \eqn{\xi_{jk}} stand for the k-th FPC score of the j-th predictor process. \code{MultiFAM} only focuses on mutiple predictor processes case. In fact, the case of univariate predictor is the same with functional additive model proposed by Mueller and Yao (2008). Especially in this development, one can designate an estimation support of additive models when the additive modeling is only allowed over restricted intervals or one is interested in the modeling over the support (see Han et al., 2016).
#'
#' @return A list containing the following fields:
#' \item{xi}{An \emph{N} by \eqn{(K_1 + ... + K_d)} matrix whose column vectors consist of FPC score grid vectors that each additive component functional is evluated.}
#' \item{SBFit}{An \emph{N} by \eqn{(K_1 + ... + K_d)} matrix whose column vectors consist of the smooth backfitting component function estimators at the given estimation points.}
#' \item{mY}{A scalar of centered part of the regression model.}
#' \item{NW}{An \emph{N} by \eqn{(K_1 + ... + K_d)} matrix whose column vectors consist of the Nadaraya-Watson marginal regression function estimators for each predictor component at the given estimation points.}
#' @examples
#' set.seed(1000)
#' 
#' library(MASS)
#' 
#' trapzRcpp <- fdapace:::trapzRcpp
#' 
#' f11 <- function(t) t
#' f12 <- function(t) 2*cos(2*pi*t/4)
#' f21 <- function(t) 1.5*sin(2*pi*t/4)
#' f22 <- function(t) 1.5*atan(2*pi*t/4)
#' 
#' n<-100
#' 
#' sig <- matrix(c(2.0, 0.0, 0.5, -.2,
#'                 0.0, 1.2, -.2, 0.3,
#'                 0.5, -.2, 1.7, 0.0,
#'                 -.2, 0.3, 0.0, 1.0),
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
#' xi <- matrix(rep(seq(-2,2,length.out=101),4),nrow=101,ncol=4)
#' sbf <- MultiFAM(Y=Y,Xi=Xi,xi=xi)
#' 
#' par(mfrow=c(2,2))
#' j <- 1
#' p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#' g11 <- f11(sort(xi[,j])) - trapzRcpp(sort(xi[,j]),f11(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
#' tmpSgn <- sign(sum(g11*sbf$SBFit[,j]))
#' plot(sort(xi[,j]),g11,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi11')
#' points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
#' 
#' j <- 2
#' p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#' g12 <- f12(sort(xi[,j])) - trapzRcpp(sort(xi[,j]),f12(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
#' tmpSgn <- sign(sum(g12*sbf$SBFit[,j]))
#' plot(sort(xi[,j]),g12,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi12')
#' points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
#' 
#' j <- 3
#' p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#' g21 <- f21(sort(xi[,j])) - trapzRcpp(sort(xi[,j]),f21(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
#' tmpSgn <- sign(sum(g21*sbf$SBFit[,j]))
#' plot(sort(xi[,j]),g21,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi21')
#' points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
#' 
#' j <- 4
#' p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#' g22 <- f22(sort(xi[,j])) - trapzRcpp(sort(xi[,j]),f22(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
#' tmpSgn <- sign(sum(g22*sbf$SBFit[,j]))
#' plot(sort(xi[,j]),g22,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi22')
#' points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
#' 
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

MultiFAM <- function(Y,Xi,K='epan',FVE=NULL,xi=NULL,h=NULL,supp=NULL){
  
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
  
  supp00 <- matrix(rep(c(-2,2),d0),ncol=2,byrow=TRUE)
  
  x0 <- c()
  
  dj <- c()
  X <- c()
  supp0 <- matrix(ncol=2)
  for (j in 1:d0) {
    tmpLy <- Xi[[j]]$Ly
    tmpLt <- Xi[[j]]$Lt
    
    tmpFPCA <- FPCA(tmpLy, tmpLt, optns = list(FVEthreshold=FVE[j]))
    
    Xij <- t(t(tmpFPCA$xiEst)/sqrt(tmpFPCA$lambda))
    dj <- length(tmpFPCA$lambda)
    
    X <- cbind(X,Xij)
    
    xj0 <- seq(supp00[j,1],supp00[j,2],length.out=N)
    x0 <- cbind(x0, matrix(rep(xj0,dj),ncol=dj))
    
    
    supp0 <- t(cbind(t(supp0),matrix(rep(supp00[j,],dj),ncol=2)))
  }
  
  supp0 <- supp0[-1,]
  
  if (is.null(supp)==TRUE) {
    supp <- supp0
  }
  
  if (ncol(xi)!=ncol(X)) {
    stop('Column length of evaluation grid matrix should have the same with K1+...+Kd')
  }
  
  if (is.null(xi)==TRUE) {
    x <- x0
  } else {
    x <- xi
    
    tmpIndex <- rep(1,nrow(x))
    for (l in 1:ncol(x)) {
      tmpIndex <- tmpIndex*dunif(x[,l],supp[l,1],supp[l,2])*(supp[l,2]-supp[l,1])
    }
    tmpIndex <- which(tmpIndex==1)
    
    x <- x[tmpIndex,]
    
  }
  
  d <- ncol(X)
  
  if (K!='epan') {
    message('Epanechnikov kernel is only supported currently. It uses Epanechnikov kernel automatically')
    K<-'epan'
  }
  if (is.null(h)==TRUE) {
    # h <- rep(0.4*n^(-1/5),d)*(supp[,2]-supp[,1])
    h <- c()
    for (j in 1:d) {
      options(warn = -1) 
      h0 <- n^(-1/5)*fdapace:::GCVLwls1D1(y=Y,t=X[,j],kernel='epan',npoly=0,nder=0,dataType='Dense')$bOpt
      options(warn = 0) 
      while (h0 > (supp[j,2]-supp[j,1])/5) {
        h0 <- (supp[j,2]-supp[j,1])/6
      }
      h[j] <- h0
    }
  } else {
    h <- rep(h,dj)
  }
  
  sbf <- SBFitting(Y,x,X,h,supp=supp)
  
  return(sbf)
  
}


