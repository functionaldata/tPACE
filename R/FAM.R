#' Functional Additive Models
#'
#' Functional additive models with a single predictor process
#'
#' @param Y An \emph{n}-dimensional vector whose elements consist of scalar responses.
#' @param Lx A list of \emph{n} vectors containing the observed values for each individual. See \code{FPCA} for detail.
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order. See \code{FPCA} for detail.
#' @param nEval The number of evaluation grid points for kernel smoothing (default is 51. If it is specified as 0, then estimated FPC scores are used for evaluation grid instead of equal grid).
#' @param bwMethod The method of bandwidth selection for kernel smoothing, a positive value for designating K-fold cross-validtaion and zero for GCV (default is 5)
#' @param supp The lower and upper limits of kernel smoothing domain for studentized FPC scores (default is [-2,2]).
#' @param optns A list of options control parameters specified by list(name=value). See \code{FPCA}.
#'
#' @details \code{FAM} fits functional additive models for a scalar response and single predictor process proposed by Mueller and Yao (2007) that \deqn{E(Y | \mathbf{X}) = \sum_{k=1}^K g_{k}(\xi_{k}),} where \eqn{\xi_{k}} stand for the k-th FPC score of the the predictor process.
#'
#' @return A list containing the following fields:
#' \item{xi}{An \emph{N} by \emph{K} matrix whose column vectors consist of \emph{N} vectors of estimation points for each component function.}
#' \item{mu}{Mean estimator of \eqn{EY}}
#' \item{fam}{A \emph{N} by \emph{K} matrix whose column vectors consist of the component function estimators at the given estimation points.}
#' \item{bw}{A \emph{K}-dimensional bandwidth vector.}
#' @examples
#' 
#'set.seed(500)
#'
#'library(MASS)
#'
#'trapzRcpp <- fdapace:::trapzRcpp
#'
#'f1 <- function(t) 0.5*t
#'f2 <- function(t) 2*cos(2*pi*t/4)
#'f3 <- function(t) 2*sin(2*pi*t/4)
#'f4 <- function(t) 2*atan(2*pi*t/4)
#'
#'n<-100
#'
#'sig <- diag(c(4.0,2.0,1.5,1.2))
#'
#'scoreX <- mvrnorm(n,mu=rep(0,4),Sigma=sig)
#'Y <- f1(scoreX[,1]) + f2(scoreX[,2]) + f3(scoreX[,3]) + f4(scoreX[,4]) + rnorm(n,0,0.5)
#'
#'phi1 <- function(t) sqrt(2)*sin(2*pi*t)
#'phi2 <- function(t) sqrt(2)*sin(4*pi*t)
#'phi3 <- function(t) sqrt(2)*cos(2*pi*t)
#'phi4 <- function(t) sqrt(2)*cos(4*pi*t)
#'
#'grid <- seq(0,1,length.out=51)
#'Lt <- Lx <- list()
#'for (i in 1:n) {
#'  Lt[[i]] <- grid
#'  Lx[[i]] <- scoreX[i,1]*phi1(grid) + scoreX[i,2]*phi2(grid) + scoreX[i,3]*phi3(grid) + scoreX[i,4]*phi4(grid)
#'}
#'
#'
#'# estimation
#'fit <- FAM(Y=Y,Lx=Lx,Lt=Lt)
#'
#'xi <- fit$xi
#'
#'par(mfrow=c(2,2))
#'j <- 1
#'g1 <- f1(sort(xi[,j])) - fdapace:::trapzRcpp(sort(xi[,j]),f1(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#'tmpSgn <- sign(sum(g1*fit$fam[,j]))
#'plot(sort(xi[,j]),g1,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi1')
#'points(sort(xi[,j]),tmpSgn*fit$fam[order(xi[,j]),j],type='l')
#'
#'j <- 2
#'g2 <- f2(sort(xi[,j])) - fdapace:::trapzRcpp(sort(xi[,j]),f2(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#'tmpSgn <- sign(sum(g2*fit$fam[,j]))
#'plot(sort(xi[,j]),g2,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi2')
#'points(sort(xi[,j]),tmpSgn*fit$fam[order(xi[,j]),j],type='l')
#'
#'j <- 3
#'g3 <- f3(sort(xi[,j])) - fdapace:::trapzRcpp(sort(xi[,j]),f3(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#'tmpSgn <- sign(sum(g3*fit$fam[,j]))
#'plot(sort(xi[,j]),g3,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi3')
#'points(sort(xi[,j]),tmpSgn*fit$fam[order(xi[,j]),j],type='l')
#'
#'j <- 4
#'g4 <- f4(sort(xi[,j])) - fdapace:::trapzRcpp(sort(xi[,j]),f4(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#'tmpSgn <- sign(sum(g4*fit$fam[,j]))
#'plot(sort(xi[,j]),g4,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi4')
#'points(sort(xi[,j]),tmpSgn*fit$fam[order(xi[,j]),j],type='l')
#'
#'
#'# fitting
#'fit <- FAM(Y=Y,Lx=Lx,Lt=Lt,nEval=0)
#'yHat <- fit$mu+apply(fit$fam,1,'sum')
#'par(mfrow=c(1,1))
#'plot(yHat,Y)
#'abline(coef=c(0,1),col=2)
#'
#'
#'# R^2
#'R2 <- 1-sum((Y-yHat)^2)/sum((Y-mean(Y))^2)
#'R2
#'
#' 
#' @references
#' \cite{Mueller, H.-G. and Yao, F. (2005), "Functional additive models", JASA, Vol.103, No.484, p.1534-1544.}
#' 
#' @export

FAM <- function(Y,Lx,Lt,nEval=51,bwMethod=5,supp=c(-2,2),optns=list()){
  
  n <- length(Y)
  
  XiStd <- c()
  estLambda <- c()

  tmpFPCA <- FPCA(Lx, Lt)
    
  XiStd <- t(t(tmpFPCA$xiEst)/sqrt(tmpFPCA$lambda))
  d <- length(tmpFPCA$lambda)
  
  estLambda <- tmpFPCA$lambda
  
  N <- c()
  xiStdGrid <- c()
  if (nEval==0) {
    N <- nrow(XiStd)
    xiStdGrid <- XiStd
  } else {
    N <- nEval
    xiStdGrid <- matrix(rep(seq(supp[1],supp[2],length.out=N),d),nrow=N,ncol=d)
  }
  
  h <- c()
  for (j in 1:d) {
    if (bwMethod>0) {
      options(warn = -1) 
      h[j] <- fdapace:::CVLwls1D(y=(Y-mean(Y)),t=XiStd[,j],kernel='epan',npoly=1,nder=0,dataType='Sparse',kFolds=bwMethod)
      options(warn = 0) 
    } else {
      options(warn = -1) 
      h[j] <- fdapace:::GCVLwls1D1(y=(Y-mean(Y)),t=XiStd[,j],kernel='epan',npoly=1,nder=0,dataType='Sparse')$bOpt
      options(warn = 0) 
    }
  } 
  
  fam <- matrix(nrow=N,ncol=d)
  for (j in 1:d) {
    xiTmp <- sort(xiStdGrid[,j])
    fitTmp <- fdapace:::Lwls1D(bw=h[j],kernel_type='epan',xin=sort(XiStd[,j]),yin=(Y[order(XiStd[,j])]-mean(Y)),xout=xiTmp,npoly=1,nder=0)
    fam[,j] <- fitTmp[match(xiStdGrid[,j],xiTmp)]
  }
  yMean <- mean(Y)
  
  xiGrid <- xiStdGrid%*%diag(sqrt(estLambda))
  
  bw <- h*sqrt(estLambda)
  
  fit <- list(xi=xiGrid, mu=yMean, fam=fam, bw=bw)
  
  return(fit)
  
}



