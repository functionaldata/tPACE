#' Take derivative of an FPCA object
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA().   
#' @param derOptns A list of options to control the derivation parameters specified by \code{list(name=value)}. See `Details'. (default = NULL)
#'
#' @details Available derivation control options are 
#' \describe{
#' \item{p}{The order of the derivatives returned (default: 0, max: 2)}
#' \item{bw}{Bandwidth for smoothing the derivatives (default: p * 0.075 * S)}
#' \item{kernelType}{Smoothing kernel choice; same available types are FPCA(). default('epan')}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' derRes <- FPCAder(res)
#' @export


FPCAder <-  function (fpcaObj, derOptns = list(p=1)) {

  derOptns <- SetDerOptions(fpcaObj,derOptns = derOptns)
  p <- derOptns[['p']]
  method <- derOptns[['method']]
  bw <- derOptns[['bw']]
  kernelType <- derOptns[['kernelType']]
  
  obsGrid <- fpcaObj$obsGrid
  workGrid <- fpcaObj$workGrid
  
  if (!class(fpcaObj) %in% 'FPCA'){
    stop("FPCAder() requires an FPCA class object as basic input")
  }

  if( ! (p %in% c(1, 2))){
    stop("'FPCAder()' is requested to use a derivative order other than p = {1, 2}!")
  } 
  
  if (p == 2) {
    warning('Second derivative is experimental only.')
  }

  muDer <- lwls1d(bw, kernelType, rep(1, length(obsGrid)), obsGrid, fpcaObj$mu, obsGrid, p, p)
  phiDer <- apply(fpcaObj$phi, 2, function(phi) lwls1d(bw, kernelType, rep(1, length(workGrid)), workGrid, phi, workGrid, p, p))

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

  fpcaObj <- append(fpcaObj, list(muDer = muDer, phiDer = phiDer, 
                                  #muDer2 = muDer2, phiDer2 = phiDer2,
                                  derOptns = derOptns))
  class(fpcaObj) <- c(class(fpcaObj), 'FPCAder')
  return(fpcaObj)
}
