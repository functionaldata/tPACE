#' Create the covariance surface plot based on the results from FPCA() or FPCder().
#'
#' This function will open a new device if not instructed otherwise.
#'
#' @param fpcaObj returned object from FPCA().
#' @param inputData A list of length 2 containing the sparse/dense (unsupported yet) observations. \code{inputData} needs to contain two fields: \code{t} for a list of time points and \code{y} for a list of observations.
#' @param ... other arguments passed into matplot for plotting options
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createPathPlot(res, subset=1:5)
#' @export

createPathPlot = function(fpcaObj, subset, k=NULL, inputData=fpcaObj$inputData, showObs=!is.null(inputData), ...){

  n <- dim(fpcaObj[['xiEst']])[1]
  
  if (!is.null(inputData)) {
    if (!all(c('t', 'y') %in% names(inputData))) {
      stop('inputData does not contain the required fields `t` and `y`')
    }
  } 
  if (showObs) {
    if (is.null(inputData)) {
      stop('Cannot show the sparse observations due to unspecified input data')
    } else {
      if (length(inputData[['t']]) != n)
        stop('length of inputData mismatches that in fpcaObj')
    }
  }
  
  if (missing(subset)) {
    subset <- seq_len(n)
  }
  # browser()
  workGrid <- fpcaObj[['workGrid']]
  fit <- fitted(fpcaObj, k=k)[subset, , drop=FALSE]
  
  # make a matrix with NAs for the sparse observations.
  maxN_i <- max(sapply(inputData[['t']][subset], length))
  obst <- sapply(inputData[['t']][subset], function(x) c(x, rep(NA, maxN_i - length(x))))
  obsy <- sapply(inputData[['y']][subset], function(x) c(x, rep(NA, maxN_i - length(x))))
  matplot(obst, obsy, type='p', ...)
  matplot(workGrid, t(fit), type='l', add=TRUE, ...)
 
  invisible()
}