#' Take derivative of an FPCA object
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA().   
#' @param derOptns A list of options to control the derivation parameters specified by \code{list(name=value)}. See `Details'. (default = NULL)
#'
#' @details Available derivation control options are 
#' \describe{
#' \item{method}{The method used for obtaining the derivatives. 'DPC': cite, 'FPC': cite}
#' \item{p}{The order of the derivatives returned (default: 0, max: 2)}
#' \item{bw}{Bandwidth for smoothing the derivatives (default: p * 0.1 * S). For 'DPC', bw * 2 is used for smoothing G^(1,1)(s,t)}
#' \item{kernelType}{Smoothing kernel choice; same available types are FPCA(). default('epan')}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' derRes <- FPCAder(res)
#' @export


FPCAder <-  function (fpcaObj, derOptns = list(p=1)) {

  derOptns <- SetDerOptions(fpcaObj,derOptns = derOptns)
  p <- derOptns[['p']]
  method <- derOptns[['method']]
  bw <- derOptns[['bw']]
  kernelType <- derOptns[['kernelType']]
  k <- derOptns[['k']]

  # TODO: truncated workGrid/obsGrid may not work
  obsGrid <- fpcaObj$obsGrid
  workGrid <- fpcaObj$workGrid
  nWorkGrid <- length(workGrid)
  gridSize <- workGrid[2] - workGrid[1]
  Lt <- fpcaObj[['inputData']][['Lt']]
  Ly <- fpcaObj[['inputData']][['Ly']]
  phi <- fpcaObj[['phi']]
  fittedCov <- fpcaObj[['fittedCov']]
  lambda <- fpcaObj[['lambda']]

  if (!class(fpcaObj) %in% 'FPCA'){
    stop("FPCAder() requires an FPCA class object as basic input")
  }

  if( ! (p %in% c(1, 2))){
    stop("The derivative order p should be in {1, 2}!")
  } 

  if (p == 2 && method == 'DPC') {
    stop('\'DPC\' method does not support p = 2')
  }

  if (p == 2) {
    warning('Second derivative is experimental only.')
  } 

  if (method == 'DPC') {
    # if (!derOptns$useTrue) {
    # Get mu'(t)
    xin <- unlist(Lt)
    yin <- unlist(Ly)
    ord <- order(xin)
    xin <- xin[ord]
    yin <- yin[ord]
    muDense <- Lwls1D(bw, kernelType, xin=xin, yin=yin, xout=obsGrid)
    mu1 <- Lwls1D(bw, kernelType, xin=xin, yin=yin, xout=workGrid, npoly=p + 1, nder=p)

    # Get raw covariance
    rcov <- BinRawCov(GetRawCov(Ly, Lt, obsGrid, muDense, 'Sparse', TRUE))

    # Use 1D smoothing on G(s, t) for G^(1,0)(s, t)
    if (is.null(derOptns[['G10_1D']]) || !derOptns[['G10_1D']]) {
      cov10 <- Lwls2DDeriv(bw, kernelType, xin=rcov$tPairs, yin=rcov$meanVals,
                           win=rcov$count, xout1=workGrid, xout2=workGrid,
                           npoly=1L, nder1=1L, nder2=0L)
    } else {
      tmpGrid <- seq(min(workGrid), max(workGrid), length=nWorkGrid - 1)
      cov10 <- ConvertSupport(tmpGrid, workGrid, apply(fpcaObj[['smoothedCov']], 2, diff) / gridSize)
    }
    cov11 <- Lwls2DDeriv(bw * 2, kernelType, xin=rcov$tPairs, yin=rcov$meanVals,
                         win=rcov$count, xout1=workGrid, xout2=workGrid,
                         npoly=2L, nder1=1L, nder2=1L)
    cov11 <- (cov11 + t(cov11)) / 2
    # } else { # use true values
      # muDense <- rep(0, length(obsGrid))
      # mu1 <- rep(0, nWorkGrid)
      # cov10 <- ConvertSupport(obsGrid, workGrid, Cov=cov10True, isCrossCov=TRUE)
      # cov11 <- ConvertSupport(obsGrid, workGrid, Cov=cov11True)
      # cov10T <- ConvertSupport(obsGrid, workGrid, Cov=cov10True, isCrossCov=TRUE)
      # cov11T <- ConvertSupport(obsGrid, workGrid, Cov=cov11True)
      # rgl::persp3d(workGrid, workGrid, cov10, xlab='s', ylab='t')
      # rgl::persp3d(workGrid, workGrid, cov10T, xlab='s', ylab='t')
      # rgl::persp3d(workGrid, workGrid, cov11)
      # rgl::persp3d(workGrid, workGrid, cov11T, xlab='s', ylab='t')
    # }
    # browser()
    eig <- eigen(cov11) 
    positiveInd <- eig[['values']] >= 0
    if (sum(positiveInd) == 0) {
      stop('Derivative surface is negative definite')
    }
    lambda1 <- eig[['values']][positiveInd] * gridSize
    FVE1 <- cumsum(lambda1) / sum(lambda1)
    # TODO: select number of derivative components
    FVEthreshold1 <- 0.9999
    k <- min(which(FVE1 >= FVEthreshold1))
    lambda1 <- lambda1[seq_len(k)]
    phi1 <- apply(eig[['vectors']][, positiveInd, drop=FALSE][, seq_len(k), drop=FALSE], 2, 
                  function(tmp) 
                    tmp / sqrt(trapzRcpp(as.numeric(workGrid),
                                         as.numeric(tmp^2))))
    # phi1 <- eig[['vectors']][, seq_len(ncol(phi))]

    fittedCov1 <- phi1 %*% diag(lambda1, k) %*% t(phi1)

    # rgl::persp3d(workGrid, workGrid, fittedCov1)
    # rgl::persp3d(workGrid, workGrid, cov11)
    # rgl::persp3d(workGrid, workGrid, cov10)

    # convert phi and fittedCov to obsGrid.
    zeta <- crossprod(cov10, phi1) * gridSize
    zetaObs <- ConvertSupport(workGrid, obsGrid, phi=zeta)
    # zetaObs1 <- ConvertSupport(workGrid, obsGrid, phi=t(cov10)) %*% phi1 * gridSize
    CovObs <- ConvertSupport(workGrid, obsGrid, Cov=fittedCov)
    phiObs <- ConvertSupport(workGrid, obsGrid, phi=phi)

    # conditional expectation
    if (!is.null(derOptns[['userSigma2']])) {
      sigma2 <- derOptns[['userSigma2']]
    } else {
      sigma2 <- ifelse(is.null(fpcaObj[['rho']]), fpcaObj[['sigma2']],
                       max(fpcaObj[['sigma2']], fpcaObj[['rho']]))
    }
    xi1 <- GetCEScores(Ly, Lt, list(verbose=FALSE), 
                       muDense, obsGrid, CovObs, 
                       lambda=rep(1, ncol(zetaObs)), zetaObs, 
                       sigma2)
    xiEst1 <- t(do.call(cbind, xi1['xiEst', ]))
    # xi <- GetCEScores(Ly, Lt, list(verbose=FALSE), 
                       # muDense, obsGrid, CovObs, 
                       # lambda=lambda, phiObs, 
                       # ifelse(is.null(fpcaObj[['rho']]),
                              # fpcaObj[['sigma2']], fpcaObj[['rho']]))
    # xiEst <- t(simplify2array(xi['xiEst', ], higher=FALSE))

    ret <- append(fpcaObj, list(muDer=mu1, phiDer=phi1, xiDer=xiEst1,
                                lambdaDer=lambda1, 
                                derOptns=derOptns))
  } else if (method == 'DPC1') {
    if (p != 1) {
      stop("For method = 'DPC1', p must equal to 1")
    }
    xin <- unlist(Lt)
    yin <- unlist(Ly)
    ord <- order(xin)
    xin <- xin[ord]
    yin <- yin[ord]
    muDense <- Lwls1D(bw, kernelType, xin=xin, yin=yin, xout=obsGrid)
    mu1 <- Lwls1D(bw, kernelType, xin=xin, yin=yin, xout=workGrid, npoly=p + 1, nder=p)

    # Get raw covariance
    rcov <- BinRawCov(GetRawCov(Ly, Lt, obsGrid, muDense, 'Sparse', TRUE))

    cov10 <- Lwls2DDeriv(bw, kernelType, xin=rcov$tPairs, yin=rcov$meanVals,
                         win=rcov$count, xout1=workGrid, xout2=workGrid,
                         npoly=1L, nder1=1L, nder2=0L)
    
    cov11 <- apply(cov10, 1, function(x) 
      Lwls1D(bw, kernelType, xin=workGrid, yin=x, xout=workGrid, npoly=2, nder=1)
    )
    cov11 <- (cov11 + t(cov11)) / 2

    eig <- eigen(cov11) 
    positiveInd <- eig[['values']] >= 0
    if (sum(positiveInd) == 0) {
      stop('Derivative surface is negative definite')
    }
    lambda1 <- eig[['values']][positiveInd] * gridSize
    FVE1 <- cumsum(lambda1) / sum(lambda1)

    # TODO: select number of derivative components
    FVEthreshold1 <- 0.9999
    k <- min(which(FVE1 >= FVEthreshold1))
    lambda1 <- lambda1[seq_len(k)]
    phi1 <- apply(eig[['vectors']][, positiveInd, drop=FALSE][, seq_len(k), drop=FALSE], 2, 
                  function(tmp) 
                    tmp / sqrt(trapzRcpp(as.numeric(workGrid),
                                         as.numeric(tmp^2))))

    fittedCov1 <- phi1 %*% diag(lambda1, k) %*% t(phi1)

    # convert phi and fittedCov to obsGrid.
    zeta <- crossprod(cov10, phi1) * gridSize
    zetaObs <- ConvertSupport(workGrid, obsGrid, phi=zeta)
    # zetaObs1 <- ConvertSupport(workGrid, obsGrid, phi=t(cov10)) %*% phi1 * gridSize
    CovObs <- ConvertSupport(workGrid, obsGrid, Cov=fittedCov)
    phiObs <- ConvertSupport(workGrid, obsGrid, phi=phi)

    # conditional expectation
    if (!is.null(derOptns[['userSigma2']])) {
      sigma2 <- derOptns[['userSigma2']]
    } else {
      sigma2 <- ifelse(is.null(fpcaObj[['rho']]), fpcaObj[['sigma2']],
                       max(fpcaObj[['sigma2']], fpcaObj[['rho']]))
    }
    xi1 <- GetCEScores(Ly, Lt, list(verbose=FALSE), 
                       muDense, obsGrid, CovObs, 
                       lambda=rep(1, ncol(zetaObs)), zetaObs, 
                       sigma2)
    xiEst1 <- t(do.call(cbind, xi1['xiEst', ]))

    ret <- append(fpcaObj, list(muDer=mu1, phiDer=phi1, xiDer=xiEst1,
                                lambdaDer=lambda1, 
                                derOptns=derOptns))
  } else if (method == 'FPC') {
    muDer <- Lwls1D(bw, kernelType, rep(1, nWorkGrid), workGrid, fpcaObj$mu, workGrid, p+0, nder= p)
    phiDer <- apply(phi, 2, function(phij) Lwls1D(bw, kernelType, rep(1, nWorkGrid), workGrid, phij, workGrid, p+0, nder= p))

    # muDer2<- fpcaObj$mu
    # phiDer2 <- fpcaObj$phi
    # for (i in seq_len(p)) {
    #    # derivative
    #    muDer2 <- getDerivative(y = muDer2, t = obsGrid, ord=p)
    #    phiDer2 <- apply(phiDer2, 2, getDerivative, t= workGrid, ord=p)
    #    # smooth
    #    muDer2 <- getSmoothCurve(t=obsGrid, 
    #                            ft= muDer2,
    #                            GCV = TRUE,
    #                            kernelType = kernelType, mult=2)
    #    phiDer2 <- apply(phiDer2, 2, function(x)
    #                      getSmoothCurve(t=workGrid, ft=x, GCV =TRUE, kernelType = kernelType, mult=1))
    # }

    ret <- append(fpcaObj, list(muDer = muDer, phiDer = phiDer, derOptns = derOptns))
  } else if (method == 'FPC1') {
    xin <- unlist(Lt)
    yin <- unlist(Ly)
    ord <- order(xin)
    xin <- xin[ord]
    yin <- yin[ord]
    muDense <- Lwls1D(bw, kernelType, xin=xin, yin=yin, xout=obsGrid)
    muDer <- Lwls1D(bw, kernelType, xin=xin, yin=yin, xout=workGrid, npoly=p + 1, nder=p)

    # Get raw covariance
    rcov <- BinRawCov(GetRawCov(Ly, Lt, obsGrid, muDense, 'Sparse', TRUE))

    if (p != 1) {
      stop("'FPC1' is available only for p=1")
    }
    cov10 <- Lwls2DDeriv(bw, kernelType, xin=rcov$tPairs, yin=rcov$meanVals,
                         win=rcov$count, xout1=workGrid, xout2=workGrid,
                         npoly=1L, nder1=1L, nder2=0L)
    phiDer <- cov10 %*% phi %*% diag(1 / lambda[seq_len(ncol(phi))]) * gridSize 

    ret <- append(fpcaObj, list(muDer = muDer, phiDer = phiDer, derOptns = derOptns))
  }

  class(ret) <- c('FPCAder', class(fpcaObj))
  return(ret)
}
