#' Functional Additive Models with Multiple Predictor Processes
#'
#' Smooth backfitting procedure for functional additive models with multiple predictor processes
#'
#' @param Y An \emph{n}-dimensional vector whose elements consist of scalar responses.
#' @param X A \emph{d}-dimensional list whose components consist of two lists of \emph{Ly} and \emph{Lt} containing observation times and functional covariate values for each predictor component, respectively. For details of \emph{Ly} and \emph{Lt}, see \code{FPCA} for detail.
#' @param ker A \code{function} object representing the base kernel to be used in the smooth backfitting algorithm (default is 'epan' which is the only option supported currently).
#' @param nEval The number of evaluation grid points for kernel smoothing (default is 51. If it is specified as 0, then estimated FPC scores in the training set are used for evaluation grid instead of equal grid).
#' @param XTest A \emph{d}-dimensional list for test set of functional predictors (default is NULL). If \code{XTest} is specified, then estimated FPC scores in the test set are used for evaluation grid.
#' @param bwMethod The method of initial bandwidth selection for kernel smoothing, a positive value for designating K-fold cross-validtaion and zero for GCV (default is 0)
#' @param alpha The shrinkage factor (positive number) for bandwidth selection. See Han et al. (2016) (default is 0.7).
#' @param supp The lower and upper limits of kernel smoothing domain for studentized FPC scores, which FPC scores are divided by the square roots of eigenvalues (default is [-2,2]).
#' @param optnsList A \emph{d}-dimensional list whose components consist of \code{optns} for each predictor component, respectively. (default is NULL which assigns the same default \code{optns} for all components as in \code{FPCA}).
#'
#' @details \code{MultiFAM} fits functional additive models for a scalar response and 
#' multiple predictor processes and implements  the smooth backfitting  algorithm provided in
#'  Han, K., Müller, H.G., Park, B.U.  (2018). Smooth backfitting for additive modeling with small errors-in-variables, 
#'  with an application to additive functional regression for multiple predictor functions. Bernoulli 24, 1233--1265.
#' 
#'   It is based on the model  \deqn{E(Y | \mathbf{X}) = \sum_{j=1}^d \sum_{k=1}^{K_j} g_{jk}(\xi_{jk}),} where \eqn{\xi_{jk}} stand for the k-th FPC scores of the j-th predictor 
#'  process. \code{MultiFAM} only is for the multiple predictor processes case.
#'  For a univariate predictor use FAM, the functional additive model (Müller and Yao 2008). 
#'  It is necessary to designate an estimation support for the additive component functions where the additive modeling is only allowed over 
#'  restricted intervals  (see Han et al., 2018).
#'
#' @return A list containing the following fields:
#' \item{mu}{A scalar for the centered regression model.}
#' \item{SBFit}{An \emph{N} by \eqn{(K_1 + ... + K_d)} matrix whose column vectors consist of the smooth backfitting component 
#' function estimators at the given \emph{N} estimation points.}
#' \item{xi}{An \emph{N} by \eqn{(K_1 + ... + K_d)} matrix whose column vectors consist of the FPC score grid vectors 
#' at which each additive component function is evaluated.}
#' \item{bw}{A \eqn{(K_1 + ... + K_d)}-dimensional bandwidth vector.}
#' \item{lambda}{A \eqn{(K_1 + ... + K_d)}-dimensional vector containing eigenvalues.}
#' \item{phi}{A \emph{d}-dimensional list whose components consist of an \emph{nWorkGrid} by \emph{K_j} matrix containing eigenfunctions, 
#' supported by \code{WorkGrid}. See \code{FPCA}.}
#' \item{workGrid}{A \emph{d}-dimensional list whose components consist of an \emph{nWorkGrid} by \emph{K_j} working grid, 
#' a regular grid on which the eigenanalysis is carried out See \code{FPCA}.}
#' @examples
#' set.seed(1000)
#' 
#' library(MASS)
#' 
#' f11 <- function(t) t
#' f12 <- function(t) 2*cos(2*pi*t/4)
#' f21 <- function(t) 1.5*sin(2*pi*t/4)
#' f22 <- function(t) 1.5*atan(2*pi*t/4)
#' 
#' n<-100
#' N<-200
#' 
#' sig <- matrix(c(2.0, 0.0, 0.5, -.2,
#'                 0.0, 1.2, -.2, 0.3,
#'                 0.5, -.2, 1.7, 0.0,
#'                 -.2, 0.3, 0.0, 1.0),
#'               nrow=4,ncol=4)
#' 
#' scoreX <- mvrnorm(n,mu=rep(0,4),Sigma=sig)
#' scoreXTest <- mvrnorm(N,mu=rep(0,4),Sigma=sig)
#' 
#' Y <- f11(scoreX[,1]) + f12(scoreX[,2]) + f21(scoreX[,3]) + f22(scoreX[,4]) + rnorm(n,0,0.5)
#' YTest <- f11(scoreXTest[,1]) + f12(scoreXTest[,2]) + 
#' f21(scoreXTest[,3]) + f22(scoreXTest[,4]) + rnorm(N,0,0.5)
#' 
#' phi11 <- function(t) sqrt(2)*sin(2*pi*t)
#' phi12 <- function(t) sqrt(2)*sin(4*pi*t)
#' phi21 <- function(t) sqrt(2)*cos(2*pi*t)
#' phi22 <- function(t) sqrt(2)*cos(4*pi*t)
#' 
#' grid <- seq(0,1,length.out=21)
#' Lt <- Lx1 <- Lx2 <- list()
#' for (i in 1:n) {
#'   Lt[[i]] <- grid
#'   Lx1[[i]] <- scoreX[i,1]*phi11(grid) + scoreX[i,2]*phi12(grid) + rnorm(1,0,0.01)
#'   Lx2[[i]] <- scoreX[i,3]*phi21(grid) + scoreX[i,4]*phi22(grid) + rnorm(1,0,0.01)
#' }
#' 
#' LtTest <- Lx1Test <- Lx2Test <- list()
#' for (i in 1:N) {
#'   LtTest[[i]] <- grid
#'   Lx1Test[[i]] <- scoreXTest[i,1]*phi11(grid) + scoreXTest[i,2]*phi12(grid) + rnorm(1,0,0.01)
#'   Lx2Test[[i]] <- scoreXTest[i,3]*phi21(grid) + scoreXTest[i,4]*phi22(grid) + rnorm(1,0,0.01)
#' }
#' 
#' X1 <- list(Ly=Lx1, Lt=Lt)
#' X2 <- list(Ly=Lx2, Lt=Lt)
#' 
#' X1Test <- list(Ly=Lx1Test, Lt=LtTest)
#' X2Test <- list(Ly=Lx2Test, Lt=LtTest)
#' 
#' X <- list(X1, X2)
#' XTest <- list(X1Test, X2Test)
#' 
#' # estimation
#' sbf <- MultiFAM(Y=Y,X=X)
#' 
#' xi <- sbf$xi
#' 
#' par(mfrow=c(2,2))
#' j <- 1
#' p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#' g11 <- f11(sort(xi[,j])) - 
#' trapzRcpp(sort(xi[,j]),f11(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
#' tmpSgn <- sign(sum(g11*sbf$SBFit[,j]))
#' plot(sort(xi[,j]),g11,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi11')
#' points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
#' legend('top',c('true','SBF'),col=c(2,1),lwd=2,bty='n',horiz=TRUE)
#' 
#' j <- 2
#' p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#' g12 <- f12(sort(xi[,j])) - 
#' trapzRcpp(sort(xi[,j]),f12(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
#' tmpSgn <- sign(sum(g12*sbf$SBFit[,j]))
#' plot(sort(xi[,j]),g12,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi12')
#' points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
#' legend('top',c('true','SBF'),col=c(2,1),lwd=2,bty='n',horiz=TRUE)
#' 
#' j <- 3
#' p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#' g21 <- f21(sort(xi[,j])) - 
#' trapzRcpp(sort(xi[,j]),f21(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
#' tmpSgn <- sign(sum(g21*sbf$SBFit[,j]))
#' plot(sort(xi[,j]),g21,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi21')
#' points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
#' legend('top',c('true','SBF'),col=c(2,1),lwd=2,bty='n',horiz=TRUE)
#' 
#' j <- 4
#' p0 <- trapzRcpp(sort(xi[,j]),dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))
#' g22 <- f22(sort(xi[,j])) - 
#' trapzRcpp(sort(xi[,j]),f22(sort(xi[,j]))*dnorm(sort(xi[,j]),0,sqrt(sig[j,j])))/p0
#' tmpSgn <- sign(sum(g22*sbf$SBFit[,j]))
#' plot(sort(xi[,j]),g22,type='l',col=2,ylim=c(-2.5,2.5),xlab='xi22')
#' points(sort(xi[,j]),tmpSgn*sbf$SBFit[order(xi[,j]),j],type='l')
#' legend('top',c('true','SBF'),col=c(2,1),lwd=2,bty='n',horiz=TRUE)
#' 
#' \donttest{
#' # fitting
#' sbf <- MultiFAM(Y=Y,X=X,nEval=0)
#' yHat <- sbf$mu+apply(sbf$SBFit,1,'sum')
#' plot(yHat,Y)
#' abline(coef=c(0,1),col=2)
#' 
#' 
#' # R^2
#' R2 <- 1-sum((Y-yHat)^2)/sum((Y-mean(Y))^2)
#' R2
#' 
#' 
#' # prediction
#' sbf <- MultiFAM(Y=Y,X=X,XTest=XTest)
#' yHat <- sbf$mu+apply(sbf$SBFit,1,'sum')
#' plot(yHat,YTest)
#' abline(coef=c(0,1),col=2)
#' }
#' @references
#' \cite{Mammen, E., Linton, O. and Nielsen, J. (1999), "The existence and asymptotic properties of a backfitting projection algorithm under weak conditions", Annals of Statistics, Vol.27, No.5, p.1443-1490.}
#'
#' \cite{Mammen, E. and Park, B. U. (2006), "A simple smooth backfitting method for additive models", Annals of Statistics, Vol.34, No.5, p.2252-2271.}
#'
#' \cite{Müller, H.-G. and Yao, F. (2008), "Functional additive models", Journal of the American Statistical Association, Vol.103, No.484, p.1534-1544.}
#'
#' \cite{Han, K., Müller, H.-G. and Park, B. U. (2016), "Smooth backfitting for additive modeling with small errors-in-variables, with an application to additive functional regression for multiple predictor functions", Bernoulli (accepted).}
#' @export

MultiFAM <- function(Y,X,ker='epan',nEval=51,XTest=NULL,bwMethod=0,alpha=0.7,supp=c(-2,2),optnsList=NULL){
  
  n <- length(Y)
  d0 <- length(X)
  
  if (is.null(optnsList)==TRUE) {
    for (j in 1:d0) {
      optnsList[[j]] <- list() 
    }
  }
  
  if (length(optnsList)==1) {
    
    for (j in 1:d0) {
      optnsList <- rep(optnsList,d0)
    }
  }
  
  N <- c()
  dj <- c()
  XiStd <- xiStdGrid <- c()
  estLambda <- c()
  estEigen <- list()
  workGrid <- list()
  for (j in 1:d0) {
    tmpLy <- X[[j]]$Ly
    tmpLt <- X[[j]]$Lt
    
    tmpFPCA <- FPCA(tmpLy, tmpLt, optns = optnsList[j])
    
    XijStd <- t(t(tmpFPCA$xiEst)/sqrt(tmpFPCA$lambda))
    dj[j] <- length(tmpFPCA$lambda)
    
    estLambda <- c(estLambda,tmpFPCA$lambda)
    
    estEigen[[j]] <- tmpFPCA$phi
    workGrid[[j]] <- tmpFPCA$workGrid
      
    XiStd <- cbind(XiStd,XijStd)
    
    if (is.null(XTest)==TRUE) {
      
      if (nEval==0) {
        N <- nrow(XiStd)
        xiStdGrid <- XiStd
      } else {
        N <- nEval
        xiStdGrid <- cbind(xiStdGrid,matrix(rep(seq(supp[1],supp[2],length.out=N),dj[j]),nrow=N,ncol=dj[j]))
      }
      
    } else {
      predobj <- predict.FPCA(tmpFPCA,newLy=XTest[[j]]$Ly,newLt=XTest[[j]]$Lt,K=dj[j])
      
      XijStdTest <- predobj$scores%*%diag(1/sqrt(tmpFPCA$lambda))
        
      N <- nrow(XijStdTest)
      xiStdGrid <- cbind(xiStdGrid,XijStdTest)
    }
    
  }
  d <- sum(dj)
  
  if (d > n) {
    stop('Too many FPC scores in the model. Try first two FPC scores by imposing optnsList=list(maxK=2)')
  }
  
  if (ker!='epan') {
    stop('Epanechnikov kernel is only supported currently. It uses Epanechnikov kernel automatically')
  } 
  
  h <- c()
  for (j in 1:d) {
    if (bwMethod>0) {
      h0 <- suppressWarnings(CVLwls1D(y=Y,t=XiStd[,j],kernel='epan',npoly=0,nder=0,dataType='Sparse',kFolds=bwMethod))
    } else {
      h0 <- suppressWarnings(GCVLwls1D1(yy=Y,tt=XiStd[,j],kernel='epan',npoly=0,nder=0,dataType='Sparse')$bOpt)
    }
    h[j] <- alpha*h0
  }
  
  supp <- matrix(rep(supp,c(d,d)),nrow=d,ncol=2)

  sbf <- SBFitting(Y,xiStdGrid,XiStd,h,K=ker,supp=supp)
  
  xiGrid <- xiStdGrid%*%diag(sqrt(estLambda))
  bw <- h*sqrt(estLambda)
  phi <- estEigen
  
  sbf <- list(mu=sbf$mY,SBFit=sbf$SBFit,xi=xiGrid,bw=bw,lambda=estLambda,phi=phi,workGrid=workGrid)
  
  return(sbf)
  
}


