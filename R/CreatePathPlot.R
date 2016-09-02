#' Create the fitted sample path plot based on the results from FPCA().
#'
#' @param fpcaObj Returned object from FPCA().
#' @param subset A vector of indices or a logical vector for subsetting the
#' observations.
#' @param K The number of components to reconstruct the fitted sample paths.
#' @param inputData A list of length 2 containing the sparse/dense
#' (unsupported yet) observations. \code{inputData} needs to contain two
#' fields: \code{Lt} for a list of time points and \code{Ly} for a list of
#' observations. Default to the `inputData` field within `fpcaObj`.
#' @param showObs Whether to plot the original observations for each subject.
#' @param obsOnly Whether to show only the original curves.
#' @param showMean Whether to plot the mean function as a bold solid curve.
#' @param derOptns A list of options to control derivation parameters; see `fitted.FPCA'. (default = NULL)
#' @param ... other arguments passed into matplot for plotting options
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan',
#'             verbose=TRUE))
#' CreatePathPlot(res, subset=1:5)
#'
#' # CreatePathPlot has a lot of usages:
#' \dontrun{
#' CreatePathPlot(res)
#' CreatePathPlot(res, 1:20)
#' CreatePathPlot(res, 1:20, showObs=FALSE)
#' CreatePathPlot(res, 1:20, showMean=TRUE, showObs=FALSE)
#' CreatePathPlot(res, 1:20, obsOnly=TRUE)
#' CreatePathPlot(res, 1:20, obsOnly=TRUE, showObs=FALSE)
#' CreatePathPlot(inputData=sampWiener, subset=1:20, obsOnly=TRUE)}
#' 
#' @export

CreatePathPlot = function(fpcaObj, subset, K=NULL,
                          inputData=fpcaObj[['inputData']],
                          showObs=!is.null(inputData), 
                          obsOnly=FALSE, showMean=FALSE,
                          derOptns = list(p=0), ...) {
  
  if (missing(fpcaObj)) {
    showFit <- FALSE
    n <- length(inputData[['Lt']])
    isDer <- FALSE
  } else {
    isDer <- 'FPCAder' %in% class(fpcaObj) || (!is.null(derOptns[['p']]) &&
                                               derOptns[['p']] >= 1)
    showFit <- !obsOnly
    n <- dim(fpcaObj[['xiEst']])[1]
  }

  inargs <- list(...)
  if (!is.null(inargs[['k']])) {
    K <- inargs[['k']]
    inargs[['k']] <- NULL
    warning("specifying 'k' is deprecated. Use 'K' instead!")
  }
  
  if (isDer && missing(showObs)) {
  # makes no sense to show original observations with derivatives.
    showObs <- FALSE 
  }

  if (!is.null(inputData)) {
    if (!all(c('Lt', 'Ly') %in% names(inputData))) {
      stop('inputData does not contain the required fields `Lt` and `Ly`')
    }
  } 

  if (showObs) {
    if (is.null(inputData)) {
      stop('Cannot show the sparse observations due to unspecified input data')
    } else {
      if (length(inputData[['Lt']]) != n)
        stop('length of inputData mismatches that in fpcaObj')
    }
  }
  
  if (missing(subset)) {
    subset <- seq_len(n)
  }

  if (!missing(fpcaObj)) {
    workGrid <- fpcaObj[['workGrid']]
  } else {
    workGrid <- NA
  }

  if (showFit) {
    fit <- fitted(fpcaObj, K=K, derOptns = derOptns)[subset, , drop=FALSE]
  }   

  defaultColPalette = palette()
  args1 <- list( xlab= 's', ylab= ' ',col = defaultColPalette, pch=1)    
  args1[names(inargs)] <- inargs
  
  #matplot(obst, obsy, type='p',...)
  #args2 = list (x = obst, y = obsy, type='p' )
  
  plotx <- ploty <- numeric(0)

  if( showObs || obsOnly ) {
    # make a matrix with NAs for the sparse observations.
    maxN_i <- max(sapply(inputData[['Lt']][subset], length))
    obst <- sapply(inputData[['Lt']][subset], function(x) c(x, rep(NA, maxN_i - length(x))))
    obsy <- sapply(inputData[['Ly']][subset], function(x) c(x, rep(NA, maxN_i - length(x))))
    plotx <- c(plotx, t(obst))
    ploty <- c(ploty, t(obsy))
  }

  if (showFit) {
    plotx <- c(plotx, rep(workGrid, nrow(fit)))
    ploty <- c(ploty, t(fit))
  }

  # Make canvas
  do.call(plot, c(list(x=plotx, y=ploty, type='n'), args1))

  if (obsOnly) {
    do.call(matplot, c(list(x=obst, y=obsy, type='l', add=TRUE), args1))
  } 
  if (showObs) {
    do.call(matplot, c(list(x=obst, y=obsy, type='p', add=TRUE), args1))
  }
  if (showFit) { # plot fitted curves
    do.call(matplot, c(list(x=workGrid, y=t(fit), type='l', add=TRUE ), args1))
  }
  if (showMean) {
    if (!isDer) {
      meanCurve <- fpcaObj[['mu']]
    } else { # isDer
      meanCurve <- fpcaObj[['muDer']]
    }
    lines(workGrid, meanCurve, lty=1, lwd=2)
  }
   
  invisible()
}
