#' Functional Principal Component Analysis
#' 
#' FPCA for dense or sparse functional data. 
#' 
#' @param y A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'
#' @details Available control options are 
#' \describe{
#' \item{bwcov}{bandwidth value for covariance function; positive numeric - default: determine automatically based on 'bwcovGcv'}
#' \item{bwcovGcv}{bandwidth choice method for covariance function; 'GMeanAndGCV','CV','GCV - default: 'GMeanAndGCV'')}
#' \item{bwmu}{bandwidth choice for mean function is using CV or GCV; positive numeric - default: determine automatically based on 'bwmuGcv'}
#' \item{bwmuGcv}{bandwidth choice method for mean function; 'GMeanAndGCV','CV','GCV - default: 'GMeanAndGCV''}
#' \item{corrPlot}{make correlation plot; logical - default: FALSE}
#' \item{corrPlotType}{which type of correlation plot to show; 'Fitted', 'Raw', 'Smoothed' - default: 'Fitted'}
#' \item{dataType}{do we have sparse or dense functional data; 'Sparse', 'Dense', 'DenseWithMV', 'p>>n' - default:  determine automatically based on 'IsRegular'}
#' \item{designPlot}{make design plot; logical - default: FALSE}
#' \item{error}{assume measurement error in the dataset; logical - default: TRUE}
#' \item{FVEthreshold}{Fraction-of-Variance-Explained threshold used during the SVD of the fitted covar. function; numeric (0,1] - default: 0.9999}
#' \item{kernel}{smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" - default: "epan" for dense data else "gauss"}
#' \item{methodCov}{ method to estimate covariance; 'PACE','RARE','CrossSectional' - automatically determined, user input ignored}
#' \item{methodMu}{ method to estimate mu; 'PACE','RARE','CrossSectional' - automatically determined, user input ignored }
#' \item{maxK}{maximum number of principal components to consider; positive integer - default: min(20, N-1), N:# of curves}
#' \item{method}{method to estimate the PC scores; 'CE', 'IN' - default: 'CE'}
#' \item{newdata}{new data points to estimate; numeric - default: NULL }
#' \item{ntest1}{number of curves used for CV when choosing bandwidth; [1,N] - default: min(30, N-1), N: # of curves}
#' \item{nRegGrid}{number of support points in each direction of covariance surface; numeric - default: 51}
#' \item{numBins}{number of bins to bin the data into; positive integer > 10, default: NULL}
#' \item{screePlot}{make scree plot; logical - default: FALSE}
#' \item{selectionMethod}{the method of choosing the number of principal components K; 'FVE','AIC','BIC': default 'FVE' - only 'FVE' avaiable now/ default 'FVE')}
#' \item{shrink}{apply shrinkage to estimates of random coefficients (dense data only); logical - default: FALSE}
#' \item{outPercent}{2-element vector in [0,1] indicating the outPercent data in the boundary - default (0,1)}
#' \item{rho}{truncation threshold for the iterative residual. 'cv': choose rho by leave-one-observation out cross-validation; 'no': use the iterative sigma2 estimate - default "cv".}
#' \item{rotationCut}{2-element vector in [0,1] indicating the percent of data truncated during sigma^2 estimation; default  (0.25, 0.75))}
#' \item{useBinnedData}{'FORCE' (Enforce the # of bins), 'AUTO' (Select the # of  bins automatically), 'OFF' (Do not bin) - default: 'AUTO'}
#' \item{useBins}{Not integrated yet: whether to bin the same observed time points when 2D smoothing; logical - default: FALSE}
#' \item{userCov}{user-defined smoothed covariance function; numerical matrix - default: NULL}
#' \item{userMu}{user-defined smoothed mean function; numerical vector - default: NULL}
#' \item{verbose}{display diagnostic messages; logical - default: FALSE}
#' }
#' @return A list containing the following fields:
#' \item{sigma2}{Variance for measure error.}
#' \item{lambda}{A vector of length \emph{K} containing eigenvalues.}
#' \item{phi}{An nWorkGrid by \emph{K} matrix containing eigenfunctions, supported on workGrid.}
#' \item{xiEst}{A \emph{n} by \emph{K} matrix containing the FPC estimates.} 
#' \item{xiVar}{A list of length \emph{n}, each entry containing the variance estimates for the FPC estimates.}
#' \item{obsGrid}{The (sorted) grid points where all observation points are pooled.}
#' \item{mu}{A vector of length nObsGrid containing the mean function estimate.}
#' \item{workGrid}{A vector of length nWorkGrid. The internal regular grid on which the eigen analysis is carried on.}
#' \item{smoothedCov}{A nWorkGrid by nWorkGrid matrix of the smoothed covariance surface.}
#' \item{fittedCov}{A nWorkGrid by nWorkGrid matrix of the fitted covariance surface, which is garanteed to be nonnegative definite.}
#' \item{optns}{A list of actually used options.}
#' \item{bwMu}{The selected (or user specified) bandwidth for smoothing the mean function.}
#' \item{bwCov}{The selected (or user specified) bandwidth for smoothing the covariance function.}
#' \item{rho}{A regularizer for the measurement error variance estimate.}
#' \item{FVE}{A vector with the percentages of the total variance explained by each FPC; at most equal to the 'FVEthreshold' used.}
#' 
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createCorrPlot(res, 'Fitted')
# TODO: add reference

FPCA = function(y, t, optns = list()){
  
  # FPCA checks the data validity for the PCA function. 
  if( CheckData(y,t) ){
    cat('FPCA has stopped.')
    return(FALSE);
  }

  # FPCA sets the options structure that are still NULL
  optns = SetOptions(y, t, optns);

  
  # FPCA checks the options validity for the PCA function. 
  numOfCurves = length(y);
  if( CheckOptions(t, optns,numOfCurves) ){
    cat('FPCA has stopped.')
    return(FALSE);
  }


  if(optns$dataType == 'Dense' || optns$dataType == 'DenseWithMV'){
    # Only applicable to Dense and Regular functional data,
    # or Dense data with Missing Values 
    
    # cross sectional mean and sample covariance for dense case
    # assume no measurement error.
    ymat = List2Mat(y)

    # Define time grids
    obsGrid = sort(unique(unlist(t)))
    regGrid = obsGrid
    workGrid = seq(min(obsGrid), max(obsGrid), length.out = optns$nRegGrid)

    # get cross sectional mean and sample cov
    smcObj = GetMeanDense(ymat, optns)
    mu = smcObj$mu
    smcObj$muDense = ConvertSupport(obsGrid, workGrid, mu = mu)
    scsObj = GetCovDense(ymat, mu, optns)
    scsObj$smoothCov = ConvertSupport(obsGrid, workGrid, Cov = scsObj$smoothCov)

  } else if(optns$dataType == 'Sparse'){
    # For Sparse case
    
    # Bin the data (potentially):
    if ( optns$useBinnedData != 'OFF'){ 
        BinnedDataset <- GetBinnedDataset(y,t,optns)
        y = BinnedDataset$newy;
        t = BinnedDataset$newt; 
    }

    # Generate basic grids:
    # obsGrid:  the unique sorted pooled time points of the sample and the new data
    # regGrid: the grid of time points for which the smoothed covariance surface assumes values
    obsGrid = sort(unique( c(unlist(t), optns$newdata)));
    regGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid);


    # Get the smoothed mean curve
    smcObj = GetSmoothedMeanCurve(y, t, obsGrid, regGrid, optns)


    # Get the smoothed covariance surface
    # mu: the smoothed mean curve evaluated at times 'obsGrid'
    mu = smcObj$mu
    scsObj = GetSmoothedCovarSurface(y, t, mu, obsGrid, regGrid, optns, optns$useBins) 
    sigma2 <- scsObj$sigma2

    # workGrid: possibly truncated version of the regGrid; truncation would occur during smoothing
    workGrid <- scsObj$outGrid
  } else {
    stop('not implemented for dataType = "p>>n" yet!')
  }

  # Get the results for the eigen-analysis
  eigObj = GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, workGrid, optns)

  # Truncated obsGrid, and observations. Empty observation due to truncation has length 0.
  truncObsGrid <- obsGrid
  if (!all(abs(optns$outPercent - c(0, 1)) < .Machine$double.eps * 2)) {
    buff <- .Machine$double.eps * max(abs(truncObsGrid)) * 3
    truncObsGrid <- truncObsGrid[truncObsGrid >= min(workGrid) - buff &
                            truncObsGrid <= max(workGrid) + buff]
    tmp <- TruncateObs(y, t, truncObsGrid)
    y <- tmp$y
    t <- tmp$t
  }

  # convert phi and fittedCov to obsGrid.
  muObs <- ConvertSupport(obsGrid, truncObsGrid, mu=mu)
  phiObs <- ConvertSupport(workGrid, truncObsGrid, phi=eigObj$phi)
  CovObs <- ConvertSupport(workGrid, truncObsGrid, Cov=eigObj$fittedCov)

  # Get scores  
  if (optns$method == 'CE') {
    if (optns$rho != 'no') {
      rho <- GetRho(y, t, optns, muObs, truncObsGrid, CovObs, eigObj$lambda, phiObs, sigma2)
      sigma2 <- rho
    }

    scoresObj <- GetCEScores(y, t, optns, muObs, truncObsGrid, CovObs, eigObj$lambda, phiObs, sigma2)

  } else if (optns$method == 'IN') {
    scoresObj <- GetINScores(ymat, t, optns, muObs, eigObj$lambda, phiObs)
  }

  # Make the return object by MakeResultFPCA
  ret <- MakeResultFPCA(optns, smcObj, muObs, scsObj, eigObj,
                        scoresObj, truncObsGrid, workGrid, rho=ifelse(optns$rho =='cv', rho, NA))

  return(ret); 
}
