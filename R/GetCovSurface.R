#' Covariance Surface
#' 
#' Covariance surface estimation for dense or sparse functional data. 
#' 
#' @param Ly A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='Dense'}).
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'
#' Available control options are 
#' \describe{
#' \item{userBwCov}{The bandwidth value for the smoothed covariance function; positive numeric - default: determine automatically based on 'methodBwCov'}
#' \item{methodBwCov}{The bandwidth choice method for the smoothed covariance function; 'GMeanAndGCV' (the geometric mean of the GCV bandwidth and the minimum bandwidth),'CV','GCV' - default: 10\% of the support}
#' \item{userBwMu}{The bandwidth value for the smoothed mean function (using 'CV' or 'GCV'); positive numeric - default: determine automatically based on 'methodBwMu'}
#' \item{methodBwMu}{The bandwidth choice method for the mean function; 'GMeanAndGCV' (the geometric mean of the GCV bandwidth and the minimum bandwidth),'CV','GCV' - default: 5\% of the support} 
#' \item{dataType}{The type of design we have (usually distinguishing between sparse or dense functional data); 'Sparse', 'Dense', 'DenseWithMV', 'p>>n' - default:  determine automatically based on 'IsRegular'}
#' \item{error}{Assume measurement error in the dataset; logical - default: TRUE}
#' \item{kernel}{Smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" - default: "gauss"; dense data are assumed noise-less so no smoothing is performed. }
#' \item{kFoldMuCov}{The number of folds to be used for mean and covariance smoothing. Default: 10}
#' \item{lean}{If TRUE the 'inputData' field in the output list is empty. Default: FALSE}
#' \item{methodMuCovEst}{The method to estimate the mean and covariance in the case of dense functional data; 'cross-sectional', 'smooth' - default: 'cross-sectional'}
#' \item{nRegGrid}{The number of support points in each direction of covariance surface; numeric - default: 51}
#' \item{numBins}{The number of bins to bin the data into; positive integer > 10, default: NULL}
#' \item{rotationCut}{The 2-element vector in [0,1] indicating the percent of data truncated during sigma^2 estimation; default  (0.25, 0.75))}
#' \item{useBinnedData}{Should the data be binned? 'FORCE' (Enforce the # of bins), 'AUTO' (Select the # of  bins automatically), 'OFF' (Do not bin) - default: 'AUTO'}
#' \item{useBinnedCov}{Whether to use the binned raw covariance for smoothing; logical - default:TRUE}
#' \item{userMu}{The user-defined smoothed mean function; list of two numerical vector 't' and 'mu' of equal size, 't' must cover the support defined 'Ly' - default: NULL}
#' \item{userSigma2}{The user-defined measurement error variance. A positive scalar. If specified then no regularization is used (rho is set to 'no', unless specified otherwise). Default to `NULL`}
#' \item{useBW1SE}{Pick the largest bandwidth such that CV-error is within one Standard Error from the minimum CV-error, relevant only if methodBwMu ='CV' and/or methodBwCov ='CV'; logical - default: FALSE}
#' }
#' @return A list containing the following fields:
#' \item{cov}{A square matrix of size nWorkGrid containing the covariance surface estimate.}
#' \item{sigma2}{A numeric estimate of the variance of measurement error.}
#' \item{workGrid}{A vector of length nWorkGrid. The internal regular grid on which the covariance surface estimation is carried out.}
#' \item{bwCov}{The selected (or user specified) bandwidth for smoothing thecovariance surface.}
#' \item{optns}{A list of actually-used options relevant to the covariance surface calculation.}
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.025)
#' sampWiener <- Wiener(n, pts)
#' mu = sin(2*pi*pts)
#' sampWiener <- Sparsify(t(t(sampWiener) + mu), pts, 10)
#' res = GetCovSurface(Ly = sampWiener$Ly, Lt = sampWiener$Lt)
#' @export
GetCovSurface = function(Ly, Lt, optns = list()){
  
  # Check the data validity for further analysis
  CheckData(Ly,Lt)
  
  # Force the data to be numeric member lists and handle NA's
  #Ly <- lapply(Ly, as.numeric) 
  #Lt <- lapply(Lt, as.numeric)
  #Lt <- lapply(Lt, signif, 14)
  #inputData <- list(Ly=Ly, Lt=Lt);
  
  inputData <- HandleNumericsAndNAN(Ly,Lt);
  Ly <- inputData$Ly;
  Lt <- inputData$Lt;
  
  # Set the options structure members that are still NULL
  optns = SetOptions(Ly, Lt, optns);
  
  # Check the options validity for the PCA function. 
  numOfCurves = length(Ly);
  CheckOptions(Lt, optns,numOfCurves)
  
  # Bin the data
  if ( optns$usergrid  == FALSE & optns$useBinnedData != 'OFF'){ 
    BinnedDataset <- GetBinnedDataset(Ly,Lt,optns)
    Ly = BinnedDataset$newy;
    Lt = BinnedDataset$newt; 
    optns[['nRegGrid']] <- min(optns[['nRegGrid']],
                               BinnedDataset[['numBins']])
    inputData$Ly <- Ly
    inputData$Lt <- Lt
  }
  
  # Generate basic grids:
  # obsGrid:  the unique sorted pooled time points of the sample and the new
  # data
  # regGrid: the grid of time points for which the smoothed covariance
  # surface assumes values
  # cutRegGrid: truncated grid specified by optns$outPercent for the cov
  # functions
  obsGrid = sort(unique( c(unlist(Lt))));
  regGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid);
  outPercent <- optns$outPercent
  buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
  rangeGrid <- range(regGrid)
  minGrid <- rangeGrid[1]
  maxGrid <- rangeGrid[2]
  cutRegGrid <- regGrid[regGrid > minGrid + diff(rangeGrid) * outPercent[1] -
                          buff & 
                          regGrid < minGrid + diff(rangeGrid) * outPercent[2] +
                          buff]
  
  ymat <- List2Mat(Ly, Lt)
  
  ## Mean function
  # If the user provided a mean function use it
  firsttsMu <- Sys.time() #First time-stamp for calculation of the mean
  userMu <- optns$userMu
  if ( is.list(userMu) && (length(userMu$mu) == length(userMu$t))){
    smcObj <- GetUserMeanCurve(optns, obsGrid, regGrid, buff)
    smcObj$muDense = ConvertSupport(obsGrid, regGrid, mu = smcObj$mu)
  } else if (optns$methodMuCovEst == 'smooth') { # smooth mean
    smcObj = GetSmoothedMeanCurve(Ly, Lt, obsGrid, regGrid, optns)
  } else if (optns$methodMuCovEst == 'cross-sectional') { # cross-sectional mean
    smcObj = GetMeanDense(ymat, obsGrid, optns)
  }
  # mu: the smoothed mean curve evaluated at times 'obsGrid'
  mu <- smcObj$mu
  lasttsMu <- Sys.time()
  
  
  firsttsCov <- Sys.time() #First time-stamp for calculation of the covariance
  ## Covariance function and sigma2
  if (!is.null(optns$userCov) && optns$methodMuCovEst != 'smooth') { 
    scsObj <- GetUserCov(optns, obsGrid, cutRegGrid, buff, ymat)
  } else if (optns$methodMuCovEst == 'smooth') {
    # smooth cov and/or sigma2
    scsObj = GetSmoothedCovarSurface(Ly, Lt, mu, obsGrid, regGrid, optns,
                                     optns$useBinnedCov) 
  } else if (optns$methodMuCovEst == 'cross-sectional') {
    scsObj = GetCovDense(ymat, mu, optns)
    if (length(obsGrid) != length(cutRegGrid) || !identical(obsGrid, cutRegGrid)) {
      scsObj$smoothCov = ConvertSupport(obsGrid, cutRegGrid, Cov =
                                          scsObj$smoothCov)
    }
    scsObj$outGrid <- cutRegGrid
  }
  sigma2 <- scsObj[['sigma2']]
  bwCov = scsObj$bwCov
  workGrid <- scsObj$outGrid
  


  # Make the return object by MakeResultFPCA
  ret <- list(cov = scsObj$smoothCov,  
              sigma2 = scsObj$sigma2,
              workGrid = workGrid,
              bwCov = bwCov,
              optns = optns)
  
  # Plot the results
  # 
  if(optns$plot == TRUE){
    plot3D::persp3D(x = workGrid, y = workGrid, z = ret$cov, col = colorRampPalette( c('blue','red') )(50))
  }
  
  return(ret)
}


