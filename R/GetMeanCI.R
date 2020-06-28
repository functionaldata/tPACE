#'@title Bootstrap pointwise confidence intervals for the mean function for densely observed data.
#'@description Note that bootstrap pointwise confidence intervals do not work for sparsely observed data.
#'@param Ly A list of n vectors containing the observed values for each individual. 
#'Missing values specified by \code{NA}s are supported for dense case \code{(dataType='dense')}.
#'@param Lt A list of n vectors containing the observation time points for each 
#'individual corresponding to each element in \code{Ly}. Each vector should be sorted in ascending order.
#'@param level A number taking values in [0,1] determing the confidence level. Default: 0.95.
#'@param R An integer holding the number of bootstrap replicates. Default: 999.
#'@param optns A list of options; see \code{\link{FPCA}} for details.
#'
#'@return A list of two elements: 
#'\item{CI}{A data frame holding three variables: \code{CIgrid} --- the time grid where the CIs are evaluated; \code{lower} and \code{upper} --- the lower and upper bounds of the CIs on \code{CIgrid}.}
#'\item{level}{The confidence level of the CIs}.
#'@examples
#'n <- 30
#'tgrid <- seq(0,1,length.out=21)
#'phi1 <- function(t) sqrt(2)*sin(2*pi*t)
#'phi2 <- function(t) sqrt(2)*sin(4*pi*t)
#'Lt <- rep(list(tgrid), n)
#'Ly <- lapply(1:n, function(i){
#'  tgrid + rnorm(1,0,2) * phi1(tgrid) + rnorm(1,0,0.5) * phi2(tgrid) + rnorm(1,0,0.01)
#'  })
#'res <- GetMeanCI(Lt = Lt, Ly = Ly, level = 0.9)
#'@export

GetMeanCI <- function (Ly, Lt, level = 0.95, R = 999, optns = list()) {
  optns = SetOptions(Ly, Lt, optns)
  
  if (length(level) > 1) {
    level = level[1]
    warning("The input level has more than 1 element; only the first one is used.")
  }
  if (level < 0 | level > 1) {
    stop("Invalid input value of level.")
  }
  if (R %% 1 != 0 | R < 0) {
    stop("R should be an positive integer.")
  }
  n <- length(Ly)
  
  if (optns$methodMuCovEst == 'smooth') {
    obsGrid = sort(unique( c(unlist(Lt))))
    regGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid)
    muMat <- t(sapply(1:R, function(b) {
      ind <- sample(x = seq_len(n), size = n, replace = TRUE)
      obsGrid = sort(unique( c(unlist(Lt[ind]))))
      regGridInd <- which(regGrid <= max(obsGrid) & regGrid >= min(obsGrid))
      res <- rep(NA, length(regGrid))
      res[regGridInd] <- GetSmoothedMeanCurve(
        y = Ly[ind], t = Lt[ind], obsGrid = obsGrid,
        regGrid = regGrid[regGridInd],
        optns = optns
      )$muDense
      res
    }))
    regGridInd <- which(!apply(muMat, 2, anyNA))
    CI <- apply(muMat[,regGridInd], 2, quantile, c((1-level)/2, 1-(1-level)/2))
    CIgrid <- regGrid[regGridInd]
  } else if (optns$methodMuCovEst == 'cross-sectional') {
    ymat <- List2Mat(Ly, Lt)
    obsGrid = sort(unique( c(unlist(Lt))))
    muMat <- t(sapply(1:R, function(b) {
      ind <- sample(x = seq_len(n), size = n, replace = TRUE)
      GetMeanDense(ymat = ymat[ind,], obsGrid = obsGrid, optns = optns)$mu
    }))
    CI <- apply(muMat, 2, quantile, c((1-level)/2, 1-(1-level)/2))
    CIgrid <- obsGrid
  }
  
  if (optns$dataType == 'Sparse') {
    warning("Bootstrap CIs for the mean function may not be computed for the entire time range.")
  }
  CI <- data.frame(CIgrid, lower = CI[1,], upper = CI[2,])
  return(list(CI = CI, level = level))
}
