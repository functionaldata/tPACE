#' Sieve estimation:
#' B-spline based estimation procedure for time-varying additive models. The VCAM function can be used to perform function-to-scalar regression.
#'
#' @param Lt An \emph{n}-dimensional list of \emph{N_i}-dimensional vectors whose elements consist of longitudinal time points for each \emph{i}-th subject.
#' @param Ly An \emph{n}-dimensional list of \emph{N_i}-dimensional vectors whose elements consist of longitudinal response observations of each \emph{i}-subject corresponding to \emph{Lt}.
#' @param X An \emph{n} by \emph{d} matrix whose row vectors consist of covariate vector of additive components for each subject.
#' @param optnAdd A list of options controls B-spline parameters for additive components, specified by list(name=value). See 'Details'.
#' @param optnVc A list of options controls B-spline parameters for varying-coefficient components, specified by list(name=value). See 'Details'.
#'
#' @details \code{VCAM} provides a simple algorithm based on B-spline basis to estimate its nonparametric additive and varying-coefficient components.
#' 
#' Available control options for \emph{optnAdd} are 
#' \describe{
#' \item{nKnot}{A \emph{d}-dimensional vector which designates the number of knots for each additive component function estimation (default=10).}
#' \item{order}{A \emph{d}-dimensional vector which designates the order of B-spline basis for each additive component function estimation (default=3).}
#' \item{grid}{A \emph{N} by \emph{d} matrix whose column vector consist of evaluation grid points for each component function estimation.}
#' }
#' and control options for \emph{optnVc} are 
#' \describe{
#' \item{nKnot}{A \emph{(d+1)}-dimensional vector which designates the number of knots for overall mean function and each varying-coefficient component function estimation (default=10).}
#' \item{order}{A \emph{(d+1)}-dimensional vector which designates the order of B-spline basis for overall mean function and each varying-coefficient component function estimation (default=3).}
#' \item{grid}{A \emph{M} by \emph{(d+1)} matrix whose column vectors consist of evaluation grid points for overall mean function and each varying-coefficient component function estimation.}
#' }
#' 
#' @return A list containing the following fields:
#' \item{Lt}{The same with input given by \emph{Lt}}
#' \item{LyHat}{Fitted values corresponding to \emph{Ly}}
#' \item{phiEst}{An \emph{N} by \emph{d} matrix whose column vectors consist of estimates for each additive component function evaluated at \emph{gridX}.}
#' \item{beta0Est}{An \emph{M}-dimensional vector for overall mean function estimates evaluated at \emph{gridT}.}
#' \item{betaEst}{An \emph{M} by \emph{d} matrix whose column vectors consist of estimates for each varying-coefficient components evaluated at \emph{gridT}.}
#' \item{gridX}{The same with input given by \emph{optnAdd$grid}}
#' \item{gridT}{The same with input  given by \emph{optnVc$grid}}
#' @examples
#' 
#' library(MASS)
#' 
#' set.seed(100)
#' 
#' n <- 100
#' d <- 2
#' 
#' Lt <- list()
#' Ly <- list()
#' 
#' m <- rep(0,2)
#' S <- matrix(c(1,0.5,0.5,1),nrow=2,ncol=2)
#' X <- pnorm(mvrnorm(n,m,S))
#' 
#' beta0 <- function(t) 1.5*sin(3*pi*(t+0.5))
#' beta1 <- function(t) 3*(1-t)^2
#' beta2 <- function(t) 4*t^3
#' 
#' phi1 <- function(x) sin(2*pi*x)
#' phi2 <- function(x) 4*x^3-1
#' 
#' for (i in 1:n) {
#'   Ni <- sample(10:20,1)
#'   
#'   Lt[[i]] <- sort(runif(Ni,0,1))
#'   Ly[[i]] <- beta0(Lt[[i]]) + 
#'      beta1(Lt[[i]])*phi1(X[i,1]) + beta2(Lt[[i]])*phi2(X[i,2]) + rnorm(Ni,0,0.1)
#'   
#' }
#' 
#' 
#' vcam <- VCAM(Lt,Ly,X)
#' 
#' op <- par(no.readonly = TRUE)
#' 
#' par(mfrow=c(1,1))
#' plot(unlist(vcam$LyHat),unlist(Ly),xlab='observed Y',ylab='fitted Y')
#' abline(coef=c(0,1),col=8)
#' 
#' par(mfrow=c(1,2))
#' plot(vcam$gridX[,1],vcam$phiEst[,1],type='l',ylim=c(-1,1),xlab='x1',ylab='phi1')
#' points(vcam$gridX[,1],phi1(vcam$gridX[,1]),col=2,type='l',lty=2,lwd=2)
#' legend('topright',c('true','est'),lwd=2,lty=c(1,2),col=c(1,2))
#' 
#' plot(vcam$gridX[,2],vcam$phiEst[,2],type='l',ylim=c(-1,3),xlab='x2',ylab='phi2')
#' points(vcam$gridX[,2],phi2(vcam$gridX[,2]),col=2,type='l',lty=2,lwd=2)
#' legend('topleft',c('true','est'),lwd=2,lty=c(1,2),col=c(1,2))
#' 
#' par(mfrow=c(1,3))
#' plot(vcam$gridT,vcam$beta0Est,type='l',xlab='t',ylab='beta0')
#' points(vcam$gridT,beta0(vcam$gridT),col=2,type='l',lty=2,lwd=2)
#' legend('topright',c('true','est'),lwd=2,lty=c(1,2),col=c(1,2))
#' 
#' plot(vcam$gridT,vcam$betaEst[,1],type='l',xlab='t',ylab='beta1')
#' points(vcam$gridT,beta1(vcam$gridT),col=2,type='l',lty=2,lwd=2)
#' legend('topright',c('true','est'),lwd=2,lty=c(1,2),col=c(1,2))
#' 
#' plot(vcam$gridT,vcam$betaEst[,2],type='l',xlab='t',ylab='beta2')
#' points(vcam$gridT,beta2(vcam$gridT),col=2,type='l',lty=2,lwd=2)
#' legend('topright',c('true','est'),lwd=2,lty=c(1,2),col=c(1,2))
#' 
#' par(op)
#' 
#' @references
#' \cite{Zhang, X. and Wang, J.-L. (2015), "Varying-coefficient additive models for functional data", Biometrika, Vol.102, No.1, p.15-32.}
#'
#' @export

VCAM <- function(Lt,Ly,X,optnAdd=list(),optnVc=list()) {
  
  n <- length(Lt)
  d <- ncol(X)
  
  if (is.null(optnAdd$order)==TRUE) {
    optnAdd$order <- 3
  }
  if (is.null(optnAdd$nKnot)==TRUE) {
    optnAdd$nKnot <- min(10,ceiling(n/d-optnAdd$order-1))
  }
  
  if (length(optnAdd$order)==1) {
    optnAdd$order <- rep(optnAdd$order,d)
  }
  if (length(optnAdd$nKnot)==1) {
    optnAdd$nKnot <- rep(optnAdd$nKnot,d)
  }
  
  if (is.null(optnVc$order)==TRUE) {
    optnVc$order <- 3
  }
  if (is.null(optnVc$nKnot)==TRUE) {
    optnVc$nKnot <- min(10,ceiling(n/(d+1)-optnVc$order-1))
  }
  
  if (length(optnVc$order)==1) {
    optnVc$order <- rep(optnVc$order,d+1)
  }
  if (length(optnVc$nKnot)==1) {
    optnVc$nKnot <- rep(optnVc$nKnot,d+1)
  }
  
  if (is.null(d)==TRUE) {
    X <- matrix(X,nrow=n,ncol=d)
  }
  
  Ni <- c()
  for (i in 1:n) {
    Ni[i] <- length(Lt[[i]])
  }
  
  maxT <- max(unlist(Lt))
  minT <- min(unlist(Lt))
  
  
  if (is.null(optnAdd$grid)==TRUE) {
    tmp <- matrix(nrow=51,ncol=d)
    for (j in 1:d) {
      tmp[,j] <- seq(min(X[,j]),max(X[,j]),length.out=51)
    }
    optnAdd$grid <- tmp
  }
  
  if (is.null(optnVc$grid)==TRUE) {
    optnVc$grid <- seq(minT,maxT,length.out=51)
  }
  
  
  # B-spline estimation for additive function
  intY <- c()
  for (i in 1:n) {
    ind <- order(Lt[[i]])
    Ti <- (Lt[[i]][ind]-minT)/(maxT-minT)
    Yi <- Ly[[i]][ind]
    
    intY[i] <- trapzRcpp(Ti,Yi)
  }
  
  BSplineX <- c()
  nBasisX <- c()
  for (j in 1:d) {
    stdXj <- (X[,j]-min(X[,j]))/(max(X[,j])-min(X[,j]))
    
    nIntKnot <- optnAdd$nKnot[j] - optnAdd$order[j] - 1
    
    tmp <- GenBSpline(stdXj,nIntKnot,optnAdd$order[j])
    nBasisX[j] <- ncol(tmp)
    
    BSplineX <- cbind(BSplineX,tmp)
    
  }
  
  if (n <= ncol(BSplineX)) {
    stop('Try smaller number of knots or lower order of B-spline basis.')
  }
  
  fitAdd <- lm(intY~BSplineX)
  
  fHat <- fitAdd$coefficients[-1]
  indNa <- is.na(fHat)
  fHat[which(indNa==TRUE)] <- 0
  
  phiHat <- matrix(nrow=n,ncol=d)
  phiHat[,1] <- BSplineX[,1:nBasisX[1]]%*%fHat[1:nBasisX[1]]
  if (d > 1) {
    for (j in 2:d) {
      phiHat[,j] <- BSplineX[,(sum(nBasisX[1:(j-1)])+1):sum(nBasisX[1:j])]%*%fHat[(sum(nBasisX[1:(j-1)])+1):sum(nBasisX[1:j])]
    }
  }
  
  for (j in 1:d) {
    phiHat[,j] <- phiHat[,j] - mean(phiHat[,j])
  }
  
  
  # B-spline estimation for varying-coefficient function
  Y <- unlist(Ly)
  stdT <- (unlist(Lt)-minT)/(maxT-minT)
  
  BSplineT <- c()
  nBasisT <- c()
  for (j in 1:(d+1)) {
    nIntKnot <- optnVc$nKnot[j] - optnVc$order[j] - 1
    
    tmp <- GenBSpline(stdT,nIntKnot,optnVc$order[j])
    nBasisT[j] <- ncol(tmp)
    
    BSplineT <- cbind(BSplineT,tmp)
  }
  
  phiHatRe <- c()
  for (j in 1:d) {
    phiHatRe <- cbind(phiHatRe,rep(phiHat[,j],Ni))
  }
  
  BSplineTX <- BSplineT[,1:nBasisT[1]]
  for (j in 1:d) {
    BSplineTX <- cbind(BSplineTX,BSplineT[,(sum(nBasisT[1:j])+1):sum(nBasisT[1:(j+1)])]*phiHatRe[,j])
  }
  
  fitVc <- lm(Y~0+BSplineTX)
  
  gHat <- fitVc$coefficients
  indNa <- is.na(gHat)
  gHat[which(indNa==TRUE)] <- 0
  
  betaHat <- matrix(nrow=sum(Ni),ncol=(d+1))
  betaHat[,1] <- BSplineT[,1:nBasisT[1]]%*%gHat[1:nBasisT[1]]
  if (d > 1) {
    for (j in 1:d) {
      betaHat[,j+1] <- BSplineT[,(sum(nBasisT[1:j])+1):sum(nBasisT[1:(j+1)])]%*%gHat[(sum(nBasisT[1:j])+1):sum(nBasisT[1:(j+1)])]
    }
  }
  
  for (j in 2:(d+1)) {
    betaHat[,j] <- betaHat[,j]/trapzRcpp(sort(stdT),betaHat[order(stdT),j])
  }
  
  
  # fitted values for responses
  yHat <- fitVc$fitted.values
  
  LyHat <- list()
  LyHat[[1]] <- yHat[1:Ni[1]]
  for (i in 2:n) {
    LyHat[[i]] <- yHat[(sum(Ni[1:(i-1)])+1):sum(Ni[1:i])]
  }
  
  
  # additive function estimation
  X0 <- optnAdd$grid
  
  BSplineX0 <- c()
  nBasisX0 <- c()
  for (j in 1:d) {
    
    nIntKnot <- optnAdd$nKnot[j] - optnAdd$order[j] - 1
    
    tmp <- GenBSpline(X0[,j],nIntKnot,optnAdd$order[j])
    nBasisX0[j] <- ncol(tmp)
    
    BSplineX0 <- cbind(BSplineX0,tmp)
    
  }
  
  phiEst <- matrix(nrow=nrow(X0),ncol=d)
  phiEst[,1] <- BSplineX0[,1:nBasisX0[1]]%*%fHat[1:nBasisX0[1]]
  if (d > 1) {
    for (j in 2:d) {
      phiEst[,j] <- BSplineX0[,(sum(nBasisX0[1:(j-1)])+1):sum(nBasisX0[1:j])]%*%fHat[(sum(nBasisX0[1:(j-1)])+1):sum(nBasisX0[1:j])]
    }
  }
  
  for (j in 1:d) {
    phiEst[,j] <- phiEst[,j] - mean(phiEst[,j])
  }
  
  par(mfrow=c(1,2))
  plot(sort(X0[,1]),phiEst[order(X0[,1]),1],type='l',ylim=c(-1,1))
  plot(sort(X0[,2]),phiEst[order(X0[,2]),2],type='l',ylim=c(-1,3))

  
  # varying-coefficient function esimation
  T0 <- optnVc$grid
  
  BSplineT0 <- c()
  nBasisT0 <- c()
  for (j in 1:(d+1)) {
    nIntKnot <- optnVc$nKnot[j] - optnVc$order[j] - 1
    
    tmp <- GenBSpline(T0,nIntKnot,optnVc$order[j])
    nBasisT0[j] <- ncol(tmp)
    
    BSplineT0 <- cbind(BSplineT0,tmp)
  }
  
  betaEst <- matrix(nrow=length(T0),ncol=(d+1))
  betaEst[,1] <- BSplineT0[,1:nBasisT0[1]]%*%gHat[1:nBasisT0[1]]
  if (d > 1) {
    for (j in 1:d) {
      betaEst[,j+1] <- BSplineT0[,(sum(nBasisT0[1:j])+1):sum(nBasisT0[1:(j+1)])]%*%gHat[(sum(nBasisT0[1:j])+1):sum(nBasisT0[1:(j+1)])]
    }
  }
  
  for (j in 2:(d+1)) {
    betaEst[,j] <- betaEst[,j]/trapzRcpp(T0,betaEst[,j])
  }
  
  par(mfrow=c(1,3))
  plot(T0,betaEst[,1],type='l')
  plot(T0,betaEst[,2],type='l')
  plot(T0,betaEst[,3],type='l')

  return(list(Lt=Lt, LyHat=LyHat, phiEst=phiEst, beta0Est=betaEst[,1], betaEst=betaEst[,-1], gridT=T0, gridX=X0))
  
}
















