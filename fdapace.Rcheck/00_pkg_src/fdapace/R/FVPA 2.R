#' Functional Variance Process Analysis for dense functional data
#' 
#' @param y A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param q A scalar defining the percentile of the pooled sample residual sample used for adjustment before taking log (default: 0.1).
#' @param optns A list of options control parameters specified by \code{list(name=value)}; by default: 'error' has to be TRUE, 'FVEthreshold' is set to 0.90. See `Details in ?FPCA'.
#'
#'
#' @return A list containing the following fields:
#' \item{sigma2}{Variance estimator of the functional variance process.} 
#' \item{fpcaObjY}{FPCA object for the original data.} 
#' \item{fpcaObjR}{FPCA object for the functional variance process associated with the original data.} 
#' 
#' @examples
#' set.seed(1)
#' n <- 25
#' pts <- seq(0, 1, by=0.01)
#' sampWiener <- Wiener(n, pts)
#' # Data have to dense for FVPA to be relevant!
#' sampWiener <- Sparsify(sampWiener, pts, 101) 
#' fvpaObj <- FVPA(sampWiener$Ly, sampWiener$Lt)
#' @references
#' \cite{Hans-Georg Müller, Ulrich Stadtmüller and Fang Yao, "Functional variance processes." Journal of the American Statistical Association 101 (2006): 1007-1018}
#' @export

FVPA = function(y, t, q= 0.1, optns = list(error=TRUE, FVEthreshold = 0.9)){ 
  
  if( (q <0) || (1 < q) ){
    warning("The value of 'q' is outside [0,1]; reseting to 0.1.")
  }
  if(is.null(optns$error)){
    stop("User provided 'optns' has to provided 'error' information.")
  }
  if(is.null(optns$FVEthreshold)){
    stop("User provided 'optns' has to provided 'FVEthreshold' information.")
  }
  if(!optns$error){
    stop("FVPA is irrelevant if no error is assumed")
  }
  if (!is.null(optns[['useBinnedData']]) && optns[['useBinnedData']] == 'FORCE') {
    stop("optns$useBinnedData cannot be 'FORCE'")
  }
  optns[['useBinnedData']] <- 'OFF'
  
  fpcaObjY <- FPCA(y, t, optns)
  
  if( fpcaObjY$optns$dataType != 'Dense' ){
    
    stop(paste0("The data has to be 'Dense' for FVPA to be relevant; the current dataType is : '", fpcaObjY$optns$dataType,"'!") )
  }
  
  yFitted <- fitted(fpcaObjY);
  rawRes = GetVarianceProcess(y, t, yFitted, workGrid = fpcaObjY$workGrid, delta = 0, logarithm = FALSE )
  delta = quantile(unlist(rawRes), q);
  rm(rawRes)
  logRes = GetVarianceProcess(y, t, yFitted, workGrid = fpcaObjY$workGrid, delta = delta, logarithm = TRUE )
  
  fpcaObjR = FPCA(logRes, t, optns);
  return( list( sigma2 = fpcaObjR$sigma2, fpcaObjY = fpcaObjY, fpcaObjR = fpcaObjR))
}

GetVarianceProcess = function(y, t, yFitted, workGrid, delta =0.0, logarithm = FALSE){
  r <- list()
  for (i in 1:nrow(yFitted)){
    tempV = delta + (  y[[i]] - approx(x = workGrid, y = yFitted[i,], xout = t[[i]])$y )^2
      if(logarithm){
        r[[i]] = log(tempV);
      } else {
        r[[i]] = tempV;
      }
  }
  return(r)
}
