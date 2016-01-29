#' Take derivative of an FPCA object
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA().   
#' @param derOptns A list of options to control the derivation parameters specified by \code{list(name=value)}. See `Details'. (default = NULL)
#'
#' @details Available derivation control options are 
#' \describe{
#' \item{p}{The order of the derivatives returned (default: 0, max: 2)}
#' \item{method}{Not used}
#' \item{GCV}{Logical specifying if GCV or CV should be used to calculated the optimal bandwidth (default: FALSE)}
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
#' derRes <- deriv(res)
#' @export


deriv.FPCA <-  function (fpcaObj, derOptns = list(p=1)) {
  
  derOptns <- SetDerOptions(derOptns)
  p <- derOptns[['p']]
  method <- derOptns[['method']]
  GCV <- derOptns[['GCV']]
  kernelType <- derOptns[['kernelType']]

  if (!class(fpcaObj) %in% 'FPCA'){
    stop("deriv.FPCA() requires an FPCA class object as basic input")
  }

  if( ! (p %in% c(1, 2))){
    stop("'deriv.FPCA()' is requested to use a derivative order other than p = {1, 2}!")
  } 

  muDer <- fpcaObj$mu
  phiDer <- fpcaObj$phi

  for (i in seq_len(p)) {
    # derivative
    muDer <- getDerivative(y = muDer, t = fpcaObj$obsGrid)
    phiDer <- apply(phiDer, 2, getDerivative, t= fpcaObj$workGrid)
    # smooth
    muDer <- getSmoothCurve(t=as.vector(fpcaObj$obsGrid), 
                            ft= muDer,
                            GCV = GCV,
                            kernelType = kernelType)
    phiDer <- apply(phiDer, 2, function(x)
                      getSmoothCurve(t=as.vector(fpcaObj$workGrid), 
                                     ft=x, 
                                     GCV = GCV, 
                                     kernelType = kernelType ))
  }


  fpcaObj <- append(fpcaObj, list(muDer = muDer, 
                                  phiDer = phiDer, 
                                  derOptns = derOptns))
  class(fpcaObj) <- c(class(fpcaObj), 'FPCAder')
  fpcaObj
}
