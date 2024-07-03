#' Functional Linear Models
#'
#' Functional linear models for scalar or functional responses and scalar and/or functional predictors.
#' 
#' @param Y Either an \emph{n}-dimensional vector whose elements consist of scalar responses, or a list which contains functional responses in the form of a list LY and the time points LT at which they are observed (i.e., list(Ly = LY,Lt = LT)).
#' @param X A list of either (1) lists which contains the observed functional predictors list Lxj and the time points list Ltj at which they are observed. It needs to be of the form \code{list(list(Ly = Lx1,Lt = Lxt1),list(Ly = Lx2,Lt = Lxt2),...)}; (2) a matrix containing one or more scalar covariates; or (3) a mix of (1) and (2), in which case the scalar covariates must come after the functional ones.
#' @param XTest A list which contains the values of functional predictors for a held-out testing set.
#' @param optnsListY A list of options control parameters for the response specified by \code{list(name=value)}. See `Details' in  \code{FPCA}.
#' @param optnsListX Either (1) A list of options control parameters for the predictors specified by \code{list(name=value)}; or (2) A list of list of options, e.g. \code{list(list(name1=value1), list(name2=value2))}. See `Details' in  \code{FPCA}.
#' @param nPerm If this argument is specified, perform a permutation test to obtain the (global) p-value for the test of regression relationship between X and Y. Recommend to set to 1000 or larger if specified.
#' 
#' @return A list of the following:
#' \item{alpha}{A length-one numeric if the response Y is scalar. Or a vector of \code{length(workGridY)} of the fitted constant alpha(t) in the linear model if Y is functional.}
#' \item{betaList}{A list of fitted beta(s) vectors, one entry per functional predictor and one entry for all scalar predictors, if Y is scalar. Each of dimension \code{length(workGridX[[j]])}.
#' 
#' Or a list of fitted beta(s,t) matrices, one per predictor, if Y is functional. Each of dimension \code{length(workGridX[[j]])} times \code{length(workGridY)}.
#' }
#' \item{R2}{The functional R2}
#' \item{pv}{Permutation p-value based on the functional R2. NA if \code{nPerm} is \code{NULL}}
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
#' \item{optnsListX}{A list of list of options actually used by the FPCA for the predictor functions}
#' \item{optnsListY}{A list of options actually used by the FPCA for the response functions}
#' \item{phiY}{A \code{length(workGridY)} by k_y the estimated eigenfunctions of Y's, where k_y is number of eigenfunctions selected for Y. NULL if Y is scalar.}
#' \item{workGridY}{A vector of working grid of the response Y's. NULL if Y is scalar}
#' \item{muY}{The mean or the mean function of the response}
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

FLM <- function(Y, X, XTest=NULL, optnsListY=NULL, optnsListX=NULL, nPerm=NULL) {

  # Internally, a scalar predictor is regarded as a constant function
  
  if (length(X) == 0) {
    stop('The predictor has length 0')
  }

  isScalar <- vapply(X, is.numeric, FALSE)
  isScalarTest <- vapply(XTest, is.numeric, FALSE)
  isFunctional  <- vapply(X, is.list, FALSE)  

  if (!all(xor(isScalar, isFunctional))) {
    stop('The input `X` does not follow the required format')
  }

  if (any(isFunctional) && 
      (any(diff(which(isFunctional)) > 1) || which(isFunctional)[1] > 1)) {
    stop('All functional predictors must be specified before the scalar covariates')
  }

  dFunctional <- sum(isFunctional)
  
  if (is.null(optnsListX) || length(optnsListX) == 0) {
    optnsListX <- lapply(seq_len(dFunctional), function(x) list())
  }
  
  if (is.null(optnsListY)) {
    optnsListY <- list()
  }
    
  if (!all(vapply(optnsListX, is.list, FALSE))) {
    optnsListX <- list(optnsListX)
  }

  if (length(optnsListX)==1 && dFunctional > 1) {
    optnsListX <- rep(optnsListX, dFunctional)
  }
  
  if (!is.null(XTest) && !identical(isScalar, isScalarTest)) {
    stop('The order of covariates appears to be different in `X` and `XTest`')
  }

  # Obtain the FPCA results for x
  xFPCA <- lapply(seq_len(dFunctional), function(j) {
    FPCA(X[[j]]$Ly, X[[j]]$Lt, optns = optnsListX[[j]])
  })

  estXiList <- lapply(xFPCA, `[[`, 'xiEst')
  estEigenX <- lapply(xFPCA, `[[`, 'phi')
  djFunc <- vapply(estXiList, ncol, 1L)
  
  if (any(isScalar)) {
    scalarX <- do.call(cbind, X[isScalar])
    dTotalScalar <- ncol(scalarX) # a scalar input field may contain multiple columns

    estXiList <- c(estXiList, list(scalarX))
    # Use the identity matrix for the eigenbasis of a scalar value
    estEigenX <- c(estEigenX, list(diag(nrow=dTotalScalar)))
  } 

  if (is.null(XTest)) {
    testXiList <- estXiList
  } else {
    testXiList <- lapply(seq_len(dFunctional), function(j) {
      predobj <- predict.FPCA(xFPCA[[j]],
                              newLy=XTest[[j]]$Ly,
                              newLt=XTest[[j]]$Lt,
                              K=djFunc[j])
      predobj$scores
    })
    if (any(isScalar)) {
      testScalarX <- do.call(cbind, XTest[isScalar])
      testXiList <- c(testXiList, list(testScalarX))
    }
  }



  # FPCA for y, if applicable
  if (is.numeric(Y)) {
    resp <- matrix(Y)
    estEigenY <- matrix(1) 
    muY <- mean(Y)
  } else if (is.list(Y)) {
    yFPCA <- FPCA(Y$Ly, Y$Lt, optns = optnsListY)
    
    resp <- yFPCA$xiEst
    dk <- length(yFPCA$lambda)
    estEigenY <- yFPCA$phi
    muY <- yFPCA[['mu']]
  }

  n <- nrow(estXiList[[1]])
  N <- nrow(testXiList[[1]])
  dAll <- vapply(estXiList, ncol, 1L)
  d <- sum(dAll)
  
  if (d > n) {
    stop('Too many FPC scores in the model. 
Try to use fewer FPC scores by imposing optnsListX=list(maxK=2) or 
optnsListX=list(FVEthreshold=0.95).')
  }
  
  BR2 <- GetBR2New(resp, estXiList)
  aVec <- BR2$aVec
  bList <- BR2$bList
  R2 <- BR2$R2
    
  if (!is.null(nPerm) && nPerm > 0) {
    R2Perm <- vapply(seq_len(nPerm), function(i) {
      ind <- sample(n)
      tmp <- GetBR2New(resp[ind, , drop=FALSE], estXiList)
      tmp$R2
    }, 0.1)
    pv <- mean(R2Perm >= R2)
  } else {
    pv <- NA
  }

  # Obtain alpha and beta functions. If aVec/beta is a matrix, then rows stand for X and columns for Y
  betaList <- lapply(seq_along(estEigenX), function(j) {
    beta <- estEigenX[[j]] %*% bList[[j]] %*% t(estEigenY)
  })
  muXBeta <- Reduce(`+`, 
                    lapply(seq_len(dFunctional), function(j) {
                     crossprod(matrix(xFPCA[[j]][['mu']]), betaList[[j]]) * 
                       (xFPCA[[j]][['workGrid']][2] - xFPCA[[j]][['workGrid']][1])
                    }),
                    init=matrix(0, 1, nrow(estEigenY)))
  alpha <- matrix(aVec, nrow=1) %*% t(estEigenY) - muXBeta 

  # if (any(isScalar)) {
  #   muZBeta <- matrix(colMeans(scalarX), nrow=1) %*% betaList[[dFunctional + 1]]
  #   alpha <- alpha - muZBeta
  # }
  # browser()

  if (is.list(Y)) {
    alpha <- alpha + muY
  }

  bMat <- do.call(rbind, bList)
  estXiMat <- do.call(cbind, estXiList)
  testXiMat <- do.call(cbind, testXiList)
  yHat <- matrix(alpha + muXBeta, n, length(alpha), byrow=TRUE) + 
    estXiMat %*% bMat %*% t(estEigenY)
  yPred <- matrix(alpha + muXBeta, N, length(alpha), byrow=TRUE) + 
    testXiMat %*% bMat %*% t(estEigenY)

  res <- list(alpha=alpha,betaList=betaList, 
              muY=muY, 
              R2=R2, pv=pv, 
              yHat=yHat,yPred=yPred,
              estXi=estXiList,testXi=testXiList)

  if (is.list(Y)) {
    res[['lambdaY']] <- yFPCA[['lambda']]
    res[['phiY']] <- yFPCA[['phi']]
    res[['workGridY']] <- yFPCA[['workGrid']]
    res[['optnsListY']] <- yFPCA$optns
  }

  if (dFunctional > 0L) {
    res[['lambdaX']] <- lapply(xFPCA, `[[`, 'lambda')
    res[['phiX']] <- lapply(xFPCA, `[[`, 'phi')
    res[['workGridX']] <- lapply(xFPCA, `[[`, 'workGrid')
    res[['optnsListX']] <- lapply(xFPCA, `[[`, 'optns')
  }

  # TODO: drop dimensions, and give names to the list outputs according to the input names

  return(res)
}


GetBR2New <- function(estEtak, estXiList) {

  estEtak <- as.matrix(estEtak)

  n <- nrow(estEtak)
  kResp <- ncol(estEtak)
  estXi <- do.call(cbind, estXiList)
  dEach <- vapply(estXiList, ncol, 1L)
  cumdEach <- cumsum(dEach)
  kXi <- ncol(estXi)

  res <- lapply(seq_len(kResp), function(k) {
    flm <- lm(estEtak[, k]~estXi)
    r2 <- summary(flm)$r.squared
    a <- as.numeric(flm$coef[1])
    b <- as.numeric(flm$coef[-1])
    list(r2=r2, a=a, b=b)
  })

  # browser()
  aVec <- vapply(res, `[[`, 0, 'a')
  bList <- lapply(seq_along(estXiList), function(j) {
    if (j == 1) {
      ind <- seq_len(dEach[j])
    } else {
      ind <- seq(cumdEach[j - 1] + 1, cumdEach[j])
    }
    b <- vapply(res, function(x) x[['b']][ind], rep(0.1, dEach[j]))
    if (dEach[j] <= 1L) {
      b <- matrix(b, nrow=dEach[j])
    }
    b
  })

  # bMat <- vapply(res, `[[`, rep(0, kXi), 'b')
  # if (kXi <= 1) {
  #   bMat <- matrix(bMat, ncol=k)
  # }
  r2k <- vapply(res, `[[`, 0, 'r2')
  totalVarEtak <- apply(estEtak, 2, var) * (n - 1)
  varExplained <- sum(totalVarEtak * r2k)

  R2 <- varExplained / sum(totalVarEtak)

  list(aVec=aVec, bList=bList, R2=R2)
}
