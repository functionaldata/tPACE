#' Functional Linear Models
#'
#' Functional linear models for scalar or functional responses and functional predictors.
#' 
#' @param Y Either an \emph{n}-dimensional vector whose elements consist of scalar responses, or a list which contains functional responses in the form of a list LY and the time points LT at which they are observed (i.e., list(Ly = LY,Lt = LT)).
#' @param X A list of lists which contains the observed functional predictors list Lxj and the time points list Ltj at which they are observed. It needs to be of the form \code{list(list(Ly = Lx1,Lt = Lxt1),list(Ly = Lx2,Lt = Lxt2),...)}
#' @param XTest A list which contains the values of functional predictors for a held-out testing set.
#' @param optnsListY A list of options control parameters for the response specified by \code{list(name=value)}. See `Details' in  \code{FPCA}.
#' @param optnsListX A list of options control parameters for the predictors specified by \code{list(name=value)}. See `Details' in  \code{FPCA}.
#' 
#' @return A list of the following:
#' \item{alpha}{A length-one numeric if the response Y is scalar. Or a vector of \code{length(workGridY)} of the fitted constant alpha(t) in the linear model if Y is functional.}
#' \item{betaList}{A list of fitted beta(s) vectors, one per predictor, if Y is scalar. Each of dimension \code{length(workGridX[[j]])}.
#' 
#' Or a list of fitted beta(s,t) matrices, one per predictor, if Y is functional. Each of dimension \code{length(workGridX[[j]])} times \code{length(workGridY)}.
#' }
#' \item{yHat}{A length n vector if Y is scalar. 
#' 
#' Or an n by \code{length(workGridY)} matrix of fitted Y's from the model if Y is functional.}
#' 
#' \item{yPred}{Same as YHat if XTest is not provided. 
#' 
#' Or a length \code{length(XTest[[1]]$Ly)} vector of predicted Y's if Y is scalar.
#' 
#' Or a \code{length(XTest[[1]]$Ly)} by \code{length(workGridY)} matrix of predicted Y's if Y is functional.}
#' 
#' 
#' \item{estXi}{A list of n by k_j matrices of estimated functional principal component scores of predictors, where k_j is the number of eigenfunctions selected for each predictor.}
#' \item{testXi}{A list of n by k_j matrices of estimated functional principal component scores of predictors in XTest, with eigenfunctions fitted only with X.}
#' \item{lambdaX}{A length sum_j k_j vector of estimated eigenvalues for predictors.}
#' \item{workGridX}{A list of vectors, each is a working grid for a predictor.}
#' \item{phiY}{A \code{length(workGridY)} by k_y the estimated eigenfunctions of Y's, where k_y is number of eigenfunctions selected for Y. NULL if Y is scalar.}
#' \item{workGridY}{A vector of working grid of the response Y's. NULL if Y is scalar}
#' @examples
#' set.seed(1000)
#' 
#' library(MASS)
#'
#' ### functional covariate
#' phi1 <- function(t,k) sqrt(2)*sin(2*pi*k*t)
#' phi2 <- function(t,k) sqrt(2)*cos(2*pi*k*t)
#'
#' lambdaX <- c(1,0.7)
#' 
#' # training set
#' n <- 50
#' Xi <- matrix(rnorm(2*n),nrow=n,ncol=2)
#' 
#' denseLt <- list(); denseLy <- list()
#' sparseLt <- list(); sparseLy <- list()
#' 
#' t0 <- seq(0,1,length.out=51)
#' for (i in 1:n) {
#'  denseLt[[i]] <- t0
#'   denseLy[[i]] <- lambdaX[1]*Xi[i,1]*phi1(t0,1) + lambdaX[2]*Xi[i,2]*phi1(t0,2)
#'   
#'   ind <- sort(sample(1:length(t0),3))
#'   sparseLt[[i]] <- t0[ind]
#'   sparseLy[[i]] <- denseLy[[i]][ind]
#' }
#' 
#' denseX <- list(Ly=denseLy,Lt=denseLt)
#' sparseX <- list(Ly=sparseLy,Lt=sparseLt)
#' 
#' denseX <- list(X=denseX)
#' sparseX <- list(X=sparseX)
#' 
#' # test set
#' N <- 30
#' 
#' XiTest <- matrix(rnorm(2*N),nrow=N,ncol=2)
#' 
#' denseLtTest <- list(); denseLyTest <- list()
#' 
#' sparseLtTest <- list(); sparseLyTest <- list()
#' 
#' t0 <- seq(0,1,length.out=51)
#' for (i in 1:N) {
#'   denseLtTest[[i]] <- t0
#'   denseLyTest[[i]] <- lambdaX[1]*XiTest[i,1]*phi1(t0,1) + lambdaX[2]*XiTest[i,2]*phi1(t0,2)
#'   
#'   ind <- sort(sample(1:length(t0),5))
#'   sparseLtTest[[i]] <- t0[ind]
#'   sparseLyTest[[i]] <- denseLyTest[[i]][ind]
#' }
#' 
#' denseXTest <- list(Ly=denseLyTest,Lt=denseLtTest)
#' sparseXTest <- list(Ly=sparseLyTest,Lt=sparseLtTest)
#' 
#' denseXTest <- list(X=denseXTest)
#' sparseXTest <- list(X=sparseXTest)
#' 
#' 
#' ### scalar response
#' beta <- c(1, -1)
#' Y <- c(Xi%*%diag(lambdaX)%*%beta) + rnorm(n,0,0.5)
#' YTest <- c(XiTest%*%diag(lambdaX)%*%beta) + rnorm(N,0,0.5)
#' 
#' ## dense
#' denseFLM <- FLM(Y=Y,X=denseX,XTest=denseXTest,optnsListX=list(FVEthreshold=0.95))
#' 
#' trueBetaList <- list()
#' trueBetaList[[1]] <- cbind(phi1(denseFLM$workGridX[[1]],1),phi1(denseFLM$workGridX[[1]],2))%*%beta
#' 
#' # coefficient function estimation error (L2-norm)
#' plot(denseFLM$workGridX[[1]],denseFLM$betaList[[1]],type='l',xlab='t',ylab=paste('beta',1,sep=''))
#' points(denseFLM$workGridX[[1]],trueBetaList[[1]],type='l',col=2)
#' 
#' denseEstErr <-
#'   sqrt(trapzRcpp(denseFLM$workGridX[[1]],(denseFLM$betaList[[1]] - trueBetaList[[1]])^2))
#' denseEstErr
#' 
#' op <- par(mfrow=c(1,2))
#' plot(denseFLM$yHat,Y,xlab='fitted Y', ylab='observed Y')
#' abline(coef=c(0,1),col=8)
#' plot(denseFLM$yPred,YTest,xlab='predicted Y', ylab='observed Y')
#' abline(coef=c(0,1),col=8)
#' par(op)
#' 
#' # prediction error
#' densePredErr <- sqrt(mean((YTest - denseFLM$yPred)^2))
#' densePredErr
#' 
#' @references
#' \cite{Yao, F., Müller, H.G., Wang, J.L. (2005). Functional linear regression analysis for longitudinal data. Annals of Statistics 33, 2873--2903.}
#' \cite{Hall, P., Horowitz, J.L. (2007). Methodology and convergence rates for functional linear regression. The Annals of Statistics, 35(1), 70--91.}
#' @export



FLM <- function(Y,X,XTest=NULL,optnsListY=NULL,optnsListX=NULL){
  
  d0 <- length(X)
  
  if (is.null(optnsListX)==TRUE) {
    for (j in 1:d0) {
      optnsListX[[j]] <- list() 
    }
  }
  
  if (length(optnsListX)==1) {
    for (j in 1:d0) {
      optnsListX <- rep(optnsListX,d0)
    }
  }
  
  n <- c()
  N <- c()
  dj <- c()
  estXi <- testXi <- c()
  estLambdaX <- c()
  estEigenX <- list()
  workGridX <- list()
  for (j in 1:d0) {
    tmpLy <- X[[j]]$Ly
    tmpLt <- X[[j]]$Lt
    
    tmpFPCA <- FPCA(tmpLy, tmpLt, optns = optnsListX)
    
    estXij <- tmpFPCA$xiEst
    dj[j] <- length(tmpFPCA$lambda)
    
    estLambdaX <- c(estLambdaX,tmpFPCA$lambda)
    
    estEigenX[[j]] <- tmpFPCA$phi
    workGridX[[j]] <- tmpFPCA$workGrid
    
    estXi <- cbind(estXi,estXij)
    
    if (is.null(XTest)==TRUE) {
      testXij <- estXij
      
      testXi <- cbind(testXi,testXij)
      
    } else {
      predobj = predict.FPCA(tmpFPCA,newLy=XTest[[j]]$Ly,newLt=XTest[[j]]$Lt,K=dj[j])
      testXij <- predobj$scores
      
      testXi <- cbind(testXi,testXij)
    }
    
  }
  n <- nrow(estXi)
  N <- nrow(testXi)
  
  d <- sum(dj)
  
  
  if (d > n) {
    stop('Too many FPC scores in the model. Try a few FPC scores by imposing optnsListX=list(maxK=2) or 
         optnsListX=list(FVEthreshold=0.95).')
  }
  
  if (class(Y)=='numeric') {
    
    flm <- lm(Y~estXi)
    
    alpha <- as.numeric(flm$coef[1])
    bVec <- as.vector(flm$coef[-1])
    
    bList <- list()
    betaList <- list()
    estXiList <- testXiList <- c()
    for (j in 1:d0) {
      if (j==1) {
        bList[[j]] <- bVec[1:dj[1]]
        
        estXiList[[j]] <- estXi[,1:dj[1]]
        testXiList[[j]] <- testXi[,1:dj[1]]
      } else {
        bList[[j]] <- bVec[(sum(dj[1:(j-1)])+1):sum(dj[1:j])]
        
        estXiList[[j]] <- estXi[,(sum(dj[1:(j-1)])+1):sum(dj[1:j])]
        testXiList[[j]] <- testXi[,(sum(dj[1:(j-1)])+1):sum(dj[1:j])]
      }
      
      if (d0==1) {
        betaList[[j]] <- as.numeric(estEigenX[[j]]%*%bList[[j]])
      } else {
        betaList[[j]] <- estEigenX[[j]]%*%bList[[j]]
      }
    }
    
    # R2 <- summary(flm)$r.sq
    
    yHat <- as.numeric(flm$fitted)
    yPred <- as.numeric(alpha + testXi%*%bVec)
    
    return(list(alpha=alpha,betaList=betaList,yHat=yHat,yPred=yPred,#R2=R2,
                estXi=estXiList,testXi=testXiList,lambdaX=estLambdaX,phiX=estEigenX,workGridX=workGridX,phiY = NULL,workGridY = NULL))
  }
  
  if (class(Y)=='list') {
    
    if (is.null(optnsListY)==TRUE) {
      optnsListY <- list()
    }
    
    tmpLy <- Y$Ly
    tmpLt <- Y$Lt
    
    tmpFPCA <- FPCA(tmpLy, tmpLt, optns = optnsListY)
    
    estEtak <- tmpFPCA$xiEst
    dk <- length(tmpFPCA$lambda)
    
    estLambdaY <- tmpFPCA$lambda
    
    estEigenY <- tmpFPCA$phi
    workGridY <- tmpFPCA$workGrid
    
    alphaVec <- c()
    bMat <- matrix(nrow=d,ncol=dk)
    for (k in 1:dk) {
      flm <- lm(estEtak[,k]~estXi)
      
      alphaVec[k] <- as.numeric(flm$coef[1])
      bMat[,k] <- as.vector(flm$coef[-1])
    }      
    
    bList <- list()
    betaList <- list()
    estXiList <- testXiList <- c()
    for (j in 1:d0) {
      if (j==1) {
        bList[[j]] <- bMat[1:dj[1],]
        
        estXiList[[j]] <- estXi[,1:dj[1]]
        testXiList[[j]] <- testXi[,1:dj[1]]
      } else {
        bList[[j]] <- bMat[(sum(dj[1:(j-1)])+1):sum(dj[1:j]),]
        
        estXiList[[j]] <- estXi[,(sum(dj[1:(j-1)])+1):sum(dj[1:j])]
        testXiList[[j]] <- testXi[,(sum(dj[1:(j-1)])+1):sum(dj[1:j])]
      }
      
#      if (d0==1) {
#        betaList[[j]] <- as.numeric(estEigenX[[j]]%*%bList[[j]]%*%t(estEigenY))
#      } else {
#        betaList[[j]] <- estEigenX[[j]]%*%bList[[j]]%*%t(estEigenY)
#      }
      betaList[[j]] <- estEigenX[[j]]%*%bList[[j]]%*%t(estEigenY)
    }
    
    alpha <- c(estEigenY%*%alphaVec)
    
    yHat <- t(matrix(rep(alpha,n),nrow=length(alpha),ncol=n)) + estXi%*%bMat%*%t(estEigenY)
    yPred <-  t(matrix(rep(alpha,n),nrow=length(alpha),ncol=N)) + testXi%*%bMat%*%t(estEigenY)
    
    
    return(list(alpha=alpha,betaList=betaList,yHat=yHat,yPred=yPred,
                estXi=estXiList,testXi=testXiList,
                lambdaX=estLambdaX,phiX=estEigenX,workGridX=workGridX,
                lambdaY=estLambdaY,phiY=estEigenY,workGridY=workGridY))
    
  } 
}
