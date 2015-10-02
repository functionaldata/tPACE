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
#' \item{bwcov}{The bandwidth value for the smoothed covariance function; positive numeric - default: determine automatically based on 'bwcovMethod'}
#' \item{bwcovMethod}{The bandwidth choice method for the smoothed covariance function; 'GMeanAndGCV','CV','GCV' - default: 'GMeanAndGCV'')}
#' \item{bwmu}{The bandwidth value for the smoothed mean function (using 'CV' or 'GCV'); positive numeric - default: determine automatically based on 'bwmuMethod'}
#' \item{bwmuMethod}{The bandwidth choice method for the mean function; 'GMeanAndGCV','CV','GCV' - default: 'GMeanAndGCV''} 
#' \item{dataType}{The type of design we have (usually distinguishing between sparse or dense functional data); 'Sparse', 'Dense', 'DenseWithMV', 'p>>n' - default:  determine automatically based on 'IsRegular'}
#' \item{diagnosticsPlot}{Make diagnostics plot (design plot, mean, scree plot and first k (<=3) eigenfunctions); logical - default: FALSE}
#' \item{error}{Assume measurement error in the dataset; logical - default: TRUE}
#' \item{fitEigenValues}{Whether also to obtain a regression fit of the eigenvalues - default: FALSE}
#' \item{FVEthreshold}{Fraction-of-Variance-Explained threshold used during the SVD of the fitted covar. function; numeric (0,1] - default: 0.9999}
#' \item{kernel}{Smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" - default: "epan" for dense data else "gauss"}
#' \item{methodCov}{The method to estimate the covariance; 'PACE','RARE','CrossSectional' - automatically determined, user input ignored}
#' \item{methodMu}{The method to estimate mu; 'PACE','RARE','CrossSectional' - automatically determined, user input ignored }
#' \item{maxK}{The maximum number of principal components to consider; positive integer - default: min(20, N-1), N:# of curves}
#' \item{methodXi}{The method to estimate the PC scores; 'CE', 'IN' - default: 'CE'}
#' \item{numCVcurves}{The number of curves used for CV when choosing bandwidth; [1,N] - default: min(30, N-1), N: # of curves}
#' \item{nRegGrid}{The number of support points in each direction of covariance surface; numeric - default: 51}
#' \item{numBins}{The number of bins to bin the data into; positive integer > 10, default: NULL}
#' \item{numComponents}{The maximum number of components to return; positive integer, default: NULL}
#' \item{selectionMethod}{The method of choosing the number of principal components K; 'FVE','AIC','BIC': default 'FVE' - only 'FVE' avaiable now/ default 'FVE')}
#' \item{shrink}{Apply shrinkage to estimates of random coefficients (dense data only); logical - default: FALSE}
#' \item{outPercent}{A 2-element vector in [0,1] indicating the outPercent data in the boundary - default (0,1)}
#' \item{rho}{The truncation threshold for the iterative residual. 'cv': choose rho by leave-one-observation out cross-validation; 'no': no regularization - default "cv".}
#' \item{rotationCut}{The 2-element vector in [0,1] indicating the percent of data truncated during sigma^2 estimation; default  (0.25, 0.75))}
#' \item{useBinnedData}{Should the data be binned? 'FORCE' (Enforce the # of bins), 'AUTO' (Select the # of  bins automatically), 'OFF' (Do not bin) - default: 'AUTO'}
#' \item{useBins}{Not integrated yet: whether to bin the same observed time points when 2D smoothing; logical - default: FALSE}
#' \item{userCov}{The user-defined smoothed covariance function; numerical matrix - default: NULL}
#' \item{userMu}{The user-defined smoothed mean function; numerical vector - default: NULL}
#' \item{verbose}{Display diagnostic messages; logical - default: FALSE}
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
#' \item{fittedCov}{A nWorkGrid by nWorkGrid matrix of the fitted covariance surface, which is guaranteed to be non-negative definite.}
#' \item{optns}{A list of actually used options.}
#' \item{bwMu}{The selected (or user specified) bandwidth for smoothing the mean function.}
#' \item{bwCov}{The selected (or user specified) bandwidth for smoothing the covariance function.}
#' \item{rho}{A regularizing scalar for the measurement error variance estimate.}
#' \item{FVE}{A vector with the percentages of the total variance explained by each FPC; at most equal to the 'FVEthreshold' used.}
#' 
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- wiener(n, pts)
#' sampWiener <- sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' createCovPlot(res, 'Fitted')
#' @references
#' \cite{Yao, Fang, Hans-Georg Mueller, and Jane-Ling Wang. "Functional data analysis for sparse longitudinal data." Journal of the American Statistical Association 100, no. 470 (2005): 577-590. (Sparse data FPCA)}
#'
#' \cite{Liu, Bitao, and Hans-Georg Mueller. "Estimating derivatives for samples of sparsely observed functions, with application to online auction dynamics." Journal of the American Statistical Association 104, no. 486 (2009): 704-717. (Sparse data FPCA)}
#'
#' \cite{Castro, P. E., W. H. Lawton, and E. A. Sylvestre. "Principal modes of variation for processes with continuous sample curves." Technometrics 28, no. 4 (1986): 329-337. (Dense data FPCA)}
#' @export

FPCA = function(y, t, optns = list()){
  
  # Check the data validity for further analysis
  if( CheckData(y,t) ){
    cat('FPCA has stopped.')
    return(FALSE);
  }
  # Force the data to be list of numeric members
  y <- lapply(y, as.numeric) 
  t <- lapply(t, as.numeric)

  # Set the options structure members that are still NULL
  optns = SetOptions(y, t, optns);

  
  # Check the options validity for the PCA function. 
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
    ymat = List2Mat(y,t)

    # Define time grids
    obsGrid = sort(unique(unlist(t)))
    regGrid = obsGrid
    workGrid = seq(min(obsGrid), max(obsGrid), length.out = optns$nRegGrid) 
    
    buff <- .Machine$double.eps * max(abs(obsGrid)) * 3 
    minGrid <- min(workGrid)
    maxGrid <- max(workGrid)
    difGrid <- maxGrid - minGrid
    workGrid <- workGrid[workGrid > minGrid + difGrid * optns$outPercent[1] - buff & 
                         workGrid < minGrid + difGrid * optns$outPercent[2] + buff]
                         
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
    obsGrid = sort(unique( c(unlist(t))));
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
  if (optns$methodXi == 'CE') {
    if (optns$rho != 'no') { 
      if( length(y) > 2048 ){
        randIndx <- sample( length(y), 2048)
        rho <- GetRho(y[randIndx], t[randIndx], optns, muObs, truncObsGrid, CovObs, eigObj$lambda, phiObs, sigma2)
      } else {
        rho <- GetRho(y, t, optns, muObs, truncObsGrid, CovObs, eigObj$lambda, phiObs, sigma2)
      }
      sigma2 <- rho
    }
    scoresObj <- GetCEScores(y, t, optns, muObs, truncObsGrid, CovObs, eigObj$lambda, phiObs, sigma2)
  } else if (optns$methodXi == 'IN') {
    ymat = List2Mat(y,t)
    scoresObj <- GetINScores(ymat, t, optns, muObs, eigObj$lambda, phiObs)
  }

  if (optns$fitEigenValues) {
    fitLambda <- FitEigenValues(scsObj$rcov, workGrid, eigObj$phi, optns$numComponents)

  } else {
    fitLambda <- NULL
  }

  # Make the return object by MakeResultFPCA
  ret <- MakeResultFPCA(optns, smcObj, muObs, scsObj, eigObj,
                        scoresObj, truncObsGrid, workGrid, rho=ifelse(optns$rho =='cv', rho, NA), fitLambda=fitLambda)
    
  # Make a quick diagnostics plot     
  if(optns$diagnosticsPlot){
    createDiagnosticsPlot(t,ret);
  }

  return(ret); 
}
