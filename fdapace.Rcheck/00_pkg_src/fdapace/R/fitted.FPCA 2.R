#' Fitted functional data from FPCA object
#' 
#' Combines the zero-meaned fitted values and the interpolated mean to get the fitted values for the trajectories 
#' or the derivatives of these trajectories.
#' Estimates are given on the work-grid, not on the observation grid. Use ConvertSupport 
#' to map the estimates to your desired domain. \code{100*(1-alpha)}-percentage coverage intervals, or 
#' bands, for trajectory estimates (not derivatives) are provided. For details consult the example.
#' 
#' @param object A object of class FPCA returned by the function FPCA().   
#' @param K The integer number of the first K components used for the representation. (default: length(fpcaObj$lambda ))
#' @param derOptns A list of options to control the derivation parameters specified by \code{list(name=value)}. See `Details'. (default = NULL)
#' @param ciOptns A list of options to control the confidence interval/band specified by \code{list(name=value)}. See `Details'. (default = NULL)
#'
#' @return If \code{alpha} is \code{NULL}, \code{p>1} or functional observations are dense, an \code{n} by \code{length(workGrid)} matrix, each row of which contains a sample. Otherwise, it returns a list which consists of the following items:
#' \item{workGrid}{An evaluation grid for fitted values.}
#' \item{fitted}{An n by length(workGrid) matrix, each row of which contains a sample.}
#' \item{cvgUpper}{An n by length(workGrid) matrix, each row of which contains the upper \code{alpha}-coverage limit}
#' \item{cvgLower}{An n by length(workGrid) matrix, each row of which contains the lower \code{alpha}-coverage limit}
#' @details Available derivation control options are 
#' \describe{
#' \item{p}{The order of the derivatives returned (default: 0, max: 2)}
#' \item{method}{The method used to produce the sample of derivatives ('FPC' (default) or 'QUO'). See Liu and Müller (2009) for more details}
#' \item{bw}{Bandwidth for smoothing the derivatives (default: p * 0.10 * S)}
#' \item{kernelType}{Smoothing kernel choice; same available types are FPCA(). default('epan')}
#' }
#' @details Available confidence interval/band control options are 
#' \describe{
#' \item{alpha}{Significant level for confidence interval/band for trajectory coverage. default=0.05 (currently only work when p=0)}
#' \item{cvgMethod}{Option for trajectory coverage method between 'interval' (pointwise coverage) and 'band' (simultaneous coverage). default='band'}
#' }
#' @param ... Additional arguments
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 5:10)
#' res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' fittedY <- fitted(res, ciOptns = list(alpha=0.05))
#' 
#' workGrid <- res$workGrid
#' cvgUpper <- fittedY$cvgUpper
#' cvgLower <- fittedY$cvgLower
#' 
#' op <- par(mfrow=c(2,3))
#' ind <- sample(1:n,6)
#' for (i in 1:6) {
#'  j <- ind[i]
#'  plot(workGrid,cvgUpper[j,],type='l',ylim=c(min(cvgLower[j,]),max(cvgUpper[j,])),col=4,lty=2,
#'    xlab='t', ylab='X(t)', main=paste(j,'-th subject',sep=''))
#'  points(workGrid,cvgLower[j,],type='l',col=4,lty=2)
#'  points(res$inputData$Lt[[j]],res$inputData$Ly[[j]])
#' }
#' par(op)
#'     
#' @references
#' \cite{Yao, F., Müller, H.-G. and Wang, J.-L. "Functional data analysis for sparse longitudinal data", Journal of the American Statistical Association, vol.100, No. 470 (2005): 577-590.}
#' 
#' \cite{Liu, Bitao, and Hans-Georg Müller. "Estimating derivatives for samples of sparsely observed functions, with application to online auction dynamics." Journal of the American Statistical Association 104, no. 486 (2009): 704-717. (Sparse data FPCA)}
#' @export

fitted.FPCA <-function (object, K = NULL, derOptns = list(p=0), ciOptns = list(alpha=NULL, cvgMethod=NULL), ...) {
  ddd <- list(...)
  if (!is.null(ddd[['k']])) {
    K <- ddd[['k']]
    warning("specifying 'k' is deprecated. Use 'K' instead!")
  }
  
  derOptns <- SetDerOptions(fpcaObject = object, derOptns)
  p <- derOptns[['p']]
  method <- derOptns[['method']]
  bw <-  derOptns[['bw']] # 
  kernelType <- derOptns[['kernelType']]

  alpha <- ciOptns[['alpha']]
  if (is.null(alpha)==FALSE) {
    if (alpha <= 0 || alpha >= 1) {
      stop("'fitted.FPCA()' is requested to use a significant level between 0 and 1.")
    }
  }
  
  cvgMethod <- ciOptns[['cvgMethod']]
  if (is.null(cvgMethod)==TRUE) {
    cvgMethod <- 'band'
  } 
  
  fpcaObj <- object
  # if (class(fpcaObj) != 'FPCA'){
    # stop("fitted.FPCA() requires an FPCA class object as basic input")
  # }

  if( is.null(K) ){
    K = length( fpcaObj$lambda )
  } else {
    if( ( round(K)>=0) && ( round(K) <= length( fpcaObj$lambda ) ) ){
      K = round(K);
    } else {
      stop("'fitted.FPCA()' is requested to use more components than it currently has available. (or 'K' is smaller than 0)")
    }
  }
 
  if( ! (p %in% c(0,1,2))){
    stop("'fitted.FPCA()' is requested to use a derivative order other than p = {0,1,2}!")
  } 

  if( p < 1 ){  
    
    ZMFV = fpcaObj$xiEst[, seq_len(K), drop = FALSE] %*% t(fpcaObj$phi[, seq_len(K), drop = FALSE]);   
    IM = fpcaObj$mu 
    
    if (is.null(alpha)==TRUE || fpcaObj$optns$dataType=='Dense') {
      return( t(apply( ZMFV, 1, function(x) x + IM))) 
    } else {
      
      bwMu <- fpcaObj$bwMu
      mu = fpcaObj$mu
      phi = fpcaObj$phi
      obsGrid = fpcaObj$obsGrid
      workGrid = fpcaObj$workGrid
      lambda = fpcaObj$lambda
      # covSmooth <- fpcaObj$smoothedCov
      
      cvgUpper <- cvgLower <- matrix(nrow=nrow(fpcaObj$xiEst), ncol=length(workGrid))
      for (i in 1:nrow(fpcaObj$xiEst)) {
        xHat <- mu + ZMFV[i,]
        
        muObs <- Lwls1D(bw = bwMu, kernelType, win = rep(1,length(workGrid)),
                        xin = workGrid, yin = mu, xout = (fpcaObj$inputData)$Lt[[i]])
        
        phiObs <- apply(phi, 2, function(phiI) Lwls1D(bw = bwMu, kernelType, win = rep(1, length(workGrid)), 
                                                      xin = workGrid, yin = phiI, xout = (fpcaObj$inputData)$Lt[[i]]))
        # HI <- diag(lambda)%*%t(phiObs)
        # 
        # indI <- match(fpcaObj$inputData$Lt[[i]],obsGrid)
        # covObsI <- covSmooth[indI,indI]
        # 
        # if (is.null(fpcaObj$sigma2)==TRUE || fpcaObj$sigma2==0) {
        #   tmpSigma2 <- max(c(0.01,abs(min(diag(covObsI)))))
        #   covObsI <- covObsI + tmpSigma2*diag(1,nrow=nrow(covObsI),ncol=ncol(covObsI))
        # } else {
        #   covObsI <- covObsI + fpcaObj$sigma2*diag(1,nrow=nrow(covObsI),ncol=ncol(covObsI))
        # }
        # 
        # omegaI <- diag(lambda) - HI%*%solve(covObsI)%*%t(HI)
        
        omegaI <- fpcaObj$xiVar[[i]]
        
        tmp <- eigen(omegaI)
        tmpA <- Re(tmp$vectors)
        tmpB <- Re(tmp$values)
        tmpB[which(tmpB<0)] <- 0
        
        if (length(tmpB)==1) {
          omegaI <- tmpA*tmpB*t(tmpA)
        } else {
          omegaI <- tmpA%*%diag(tmpB)%*%t(tmpA)
        }
        
        
        if (cvgMethod=='interval') {
          cvgUpper[i,] <- xHat + stats::qnorm(1-alpha/2)*sqrt(diag(phi%*%omegaI%*%t(phi)))
          cvgLower[i,] <- xHat + stats::qnorm(alpha/2)*sqrt(diag(phi%*%omegaI%*%t(phi)))
        } else {
          cvgUpper[i,] <- xHat + sqrt(stats::qchisq(1-alpha,K)*diag(phi%*%omegaI%*%t(phi)))
          cvgLower[i,] <- xHat - sqrt(stats::qchisq(1-alpha,K)*diag(phi%*%omegaI%*%t(phi)))
        }
      }
      
      return(list(
        workGrid = workGrid,
        fitted = t(apply( ZMFV, 1, function(x) x + IM)),
        cvgUpper = cvgUpper,
        cvgLower = cvgLower
        )
      )
      
    }
  } else { #Derivative is not zero
 
    if( K > SelectK( fpcaObj, FVEthreshold=0.95, criterion='FVE')$K ){
    warning("Potentially you use too many components to estimate derivatives. \n  Consider using SelectK() to find a more informed estimate for 'K'.");
    }

    if( is.null(method) ){
      method = 'FPC'
    }

    mu = fpcaObj$mu
    phi = fpcaObj$phi
    obsGrid = fpcaObj$obsGrid
    workGrid = fpcaObj$workGrid

    if ( method == 'FPC'){
      phi = apply(phi, 2, function(phiI) Lwls1D(bw = bw, kernelType, win = rep(1, length(workGrid)), 
                                                  xin = workGrid, yin = phiI, xout = workGrid, npoly = p, nder = p))
      mu = Lwls1D(bw = bw, kernelType, win = rep(1, length(workGrid)), xin = workGrid, yin = mu, xout = workGrid, npoly = p, nder = p)
      ZMFV = fpcaObj$xiEst[, seq_len(K), drop = FALSE] %*% t(phi[, seq_len(K), drop = FALSE]);
      IM = mu 
      return( t(apply( ZMFV, 1, function(x) x + IM) ))
    }

    if( method == 'QUO'){
      impSample <- fitted(fpcaObj, K = K); # Look ma! I do recursion!
      return( t(apply(impSample, 1, function(curve) Lwls1D(bw = bw, kernelType, win = rep(1, length(workGrid)), 
                                                         xin = workGrid, yin = curve, xout = workGrid, npoly = p, nder = p))))
    } else if (method == 'DPC') {
      if (K > ncol(fpcaObj[['xiDer']])) {
        stop('fpcaObj does not contain K columns!')
      }
      return(tcrossprod(fpcaObj[['xiDer']][, seq_len(K), drop=FALSE], 
                        fpcaObj[['phiDer']][, seq_len(K), drop=FALSE]))
    }else {
      stop('You asked for a derivation scheme that is not implemented.')
    }
  }

}

getEnlargedGrid <- function(x){
  N <- length(x)
  return (  c( x[1] - 0.1 * diff(x[1:2]), x, x[N] + 0.1 * diff(x[(N-1):N])) )
}

getDerivative <- function(y, t, ord=1){  # Consider using the smoother to get the derivatives
  if( length(y) != length(t) ){
    stop("getDerivative y/t lengths are unequal.")
  }
  newt = getEnlargedGrid(t) # This is a trick to get first derivatives everywhere
  newy = Hmisc::approxExtrap(x=t, y=y, xout= newt)$y

  if (ord == 1) {
    der <- numDeriv::grad( stats::splinefun(newt, newy) , x = t )
  } else if (ord == 2) {
    der <- sapply(t, function(t0) 
                  numDeriv::hessian( stats::splinefun(newt, newy) , x = t0 )
                  )
  }

  return(der)
}

getSmoothCurve <- function(t, ft, GCV = FALSE, kernelType = 'epan', mult = 1){
  myBw = ifelse( GCV, GCVLwls1D1( yy= ft, tt =t, npoly=1, nder=0, dataType='Sparse', kernel=kernelType)[['bOpt']] ,
                      CVLwls1D(   y= ft, t = t, npoly=1, nder=0, dataType='Sparse', kernel=kernelType, kFolds = 10))
  myBw <- myBw * mult
  smoothCurve = Lwls1D(bw = myBw, kernel_type= kernelType, win = rep(1, length(t)), yin = ft, xout = t, xin= t)
  return(smoothCurve)
}




