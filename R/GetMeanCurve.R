#' Mean Curve
#' 
#' Mean curve calculation for dense or sparse functional data. 
#' 
#' @param Ly A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='Dense'}).
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'
#' Available control options are 
#' \describe{
#' \item{userBwMu}{The bandwidth value for the smoothed mean function (using 'CV' or 'GCV'); positive numeric - default: determine automatically based on 'methodBwMu'}
#' \item{methodBwMu}{The bandwidth choice method for the mean function; 'GMeanAndGCV' (the geometric mean of the GCV bandwidth and the minimum bandwidth),'CV','GCV' - default: 5\% of the support} 
#' \item{dataType}{The type of design we have (usually distinguishing between sparse or dense functional data); 'Sparse', 'Dense', 'DenseWithMV', 'p>>n' - default:  determine automatically based on 'IsRegular'}
#' \item{plot}{Plot mean curve; logical - default: FALSE}
#' \item{kernel}{Smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" - default: "gauss"; dense data are assumed noise-less so no smoothing is performed. }
#' \item{kFoldMuCov}{The number of folds to be used for mean and covariance smoothing. Default: 10}
#' \item{methodMuCovEst}{The method to estimate the mean and covariance in the case of dense functional data; 'cross-sectional', 'smooth' - default: 'cross-sectional'}
#' \item{numBins}{The number of bins to bin the data into; positive integer > 10, default: NULL}
#' \item{useBinnedData}{Should the data be binned? 'FORCE' (Enforce the # of bins), 'AUTO' (Select the # of  bins automatically), 'OFF' (Do not bin) - default: 'AUTO'}
#' \item{userMu}{The user-defined smoothed mean function; list of two numerical vector 't' and 'mu' of equal size, 't' must cover the support defined 'Ly' - default: NULL}
#' \item{useBW1SE}{Pick the largest bandwidth such that CV-error is within one Standard Error from the minimum CV-error, relevant only if methodBwMu ='CV' and/or methodBwCov ='CV'; logical - default: FALSE}
#' }
#' @return A list containing the following fields:
#' \item{mu}{A vector of length nWorkGrid containing the mean function estimate.}
#' \item{workGrid}{A vector of length nWorkGrid. The internal regular grid on which the mean estimation is carried out.}
#' \item{bwMu}{The selected (or user specified) bandwidth for smoothing the mean function.}
#' \item{optns}{A list of actually-used options relevant to the mean function calculation.}
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.025)
#' sampWiener <- Wiener(n, pts)
#' mu = sin(2*pi*pts)
#' sampWiener <- Sparsify(t(t(sampWiener) + mu), pts, 10)
#' res = GetMeanCurve(Ly = sampWiener$Ly, Lt = sampWiener$Lt, optns = list(plot = TRUE))
#' @export

GetMeanCurve = function(Ly, Lt, optns = list()){
  
  firsttsMean <- Sys.time() #First time-stamp for mean
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
  # convert mu to truncated workGrid
  workGrid <- cutRegGrid
  mu = spline(x= obsGrid, y= smcObj$mu, xout=  workGrid)$y
  
  lasttsMu <- Sys.time()
  bwMu = smcObj$bwMu
  

  # Make the return object by MakeResultFPCA
  ret <- list(mu = mu,
              workGrid = workGrid,
              bwMu = bwMu,
              optns = optns
  )

  # Plot the results
    if(optns$plot){
      plot(workGrid, mu, type='l', xlab='s',ylab='', main='Mean Function', panel.first = grid(), axes = TRUE)   
    }
  
  return(ret) 
}

