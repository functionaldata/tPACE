#' Fitted functional data for derivatives from the FPCAder object
#' 
#' Combines the zero-meaned fitted values and the mean derivative to get the fitted values for the derivative trajectories.
#' Estimates are given on the work-grid, not on the observation grid. Use ConvertSupport to map the 
#' estimates to your desired domain.
#' 
#' @param object A object of class FPCA returned by the function FPCA().   
#' @param K The integer number of the first K components used for the representation. (default: length(derObj$lambda ))
#' @param ... Additional arguments
#'
#' @return An \code{n} by \code{length(workGrid)} matrix, each row of which contains a sample.
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
# #' res <- FPCA(sampWiener$Ly, sampWiener$Lt, 
# #'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
# #' fittedY <- fitted(res)
#' @references
#' \cite{Liu, Bitao, and Hans-Georg Müller. "Estimating derivatives for samples of sparsely observed functions, with application to online auction dynamics." Journal of the American Statistical Association 104, no. 486 (2009): 704-717. (Sparse data FPCA)}
#' @export


fitted.FPCAder <-  function (object, K = NULL, ...) {
  ddd <- list(...)
  if (!is.null(ddd[['k']])) {
    K <- ddd[['k']]
    warning("specifying 'k' is deprecated. Use 'K' instead!")
  }
  
  derObj <- object
  # if (class(derObj) != 'FPCA'){
    # stop("fitted.FPCA() requires an FPCA class object as basic input")
  # }
  method <- derObj[['derOptns']][['method']]

  if (substr(method, 1, 3) == 'DPC') {
    maxK <- length(derObj[['lambdaDer']])
  } else {
    maxK <- length(derObj[['lambda']])
  }

  if( is.null(K)) {
    K <- maxK
  } else if (abs(K - round(K)) > 1e-5 || K <= 0) {
    stop("'K' needs to be a positive integer")
  } else if (K > maxK) {
      stop("'fitted.FPCAder()' is requested to use more components than it currently has available.")
  }
 
  if (substr(method, 1, 3) == 'DPC') {
    fit <- tcrossprod(derObj[['xiDer']][, seq_len(K), drop=FALSE],
                      derObj[['phiDer']][, seq_len(K), drop=FALSE])
  } else if (substr(method, 1, 3) == 'FPC') {
    fit <- tcrossprod(derObj[['xiEst']][, seq_len(K), drop=FALSE],
                      derObj[['phiDer']][, seq_len(K), drop=FALSE])
  }
  fit <- fit + matrix(derObj[['muDer']], nrow(fit), ncol(fit), byrow=TRUE)

  return(fit) 
}
