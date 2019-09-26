#' Functional Cross Covariance between longitudinal variable Y and longitudinal variable X
#' 
#' Calculate the raw and the smoothed cross-covariance between functional predictors using bandwidth bw or estimate that bw using GCV. 
#' 
#' @param Ly1 List of N vectors with amplitude information (Y)
#' @param Lt1 List of N vectors with timing information (Y)
#' @param Ymu1 Vector Q-1 Vector of length nObsGrid containing the mean function estimate (Y)
#' @param bw1 Scalar bandwidth for smoothing the cross-covariance function (if NULL it will be automatically estimated) (Y)
#' @param Ly2 List of N vectors with amplitude information (X)
#' @param Lt2 List of N vectors with timing information (X)
#' @param Ymu2 Vector Q-1 Vector of length nObsGrid containing the mean function estimate (X)
#' @param bw2 Scalar bandwidth for smoothing the cross-covariance function (if NULL it will be automatically estimated) (X)
#' @param useGAM Indicator to use gam smoothing instead of local-linear smoothing (semi-parametric option) (default: FALSE)
#' @param rmDiag Indicator to remove the diagonal element when smoothing (default: FALSE)
#' @param kern String specifying the kernel type (default: FALSE;  see ?FPCA for details)
#' @param bwRoutine String specifying the routine used to find the optimal bandwidth 'grid-search', 'bobyqa', 'l-bfgs-b' (default: 'l-bfgs-b')
#' If the variables Ly1 and Ly2 are in matrix form the data are assumed dense
#' and only the raw cross-covariance is returned. One can obtain Ymu1 and Ymu2
#' from \code{FPCA} and \code{ConvertSupport}.
#' @return A list containing:
#' \item{smoothedCC}{The smoothed cross-covariance as a matrix (currently only 51 by 51)}
#' \item{rawCC}{The raw cross-covariance as a list}
#' \item{bw}{The bandwidth used for smoothing as a vector of length 2}
#' \item{score}{The GCV score associated with the scalar used}
#' \item{smoothGrid}{The grid over which the smoothed cross-covariance is evaluated}
#' @examples
#' Ly1= list( rep(2.1,7), rep(2.1,3),2.1 );
#' Lt1 = list(1:7,1:3, 1);
#' Ly2 = list( rep(1.1,7), rep(1.1,3),1.1); 
#' Lt2 = list(1:7,1:3, 1);
#' Ymu1 = rep(55,7);
#' Ymu2 = rep(1.1,7);
#' AA<-GetCrCovYX(Ly1 = Ly1, Ly2= Ly2, Lt1=Lt1, Lt2=Lt2, Ymu1=Ymu1, Ymu2=Ymu2)
#'   
#' @references
#' \cite{Yang, Wenjing, Hans-Georg Müller, and Ulrich Stadtmüller. "Functional singular component analysis." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 73.3 (2011): 303-324}
#' @export

GetCrCovYX <- function(bw1 = NULL, bw2 = NULL, Ly1, Lt1 = NULL, Ymu1 = NULL, Ly2, Lt2 = NULL, Ymu2 = NULL, 
                       useGAM = FALSE, rmDiag=FALSE, kern='gauss', bwRoutine = 'l-bfgs-b') {
  
  if (kern != 'gauss' && (is.null(bw1) || is.null(bw2))) {
    stop('Cannot select bandwidth for non-Gaussian kernel')
  }
  
  if( !(bwRoutine %in% c('grid-search', 'bobyqa', 'l-bfgs-b') ) ){ 
    stop("'bwRoutine' argument is invalid.")
  }
  
  # If only Ly1 and Ly2 are available assume DENSE data
  if( is.matrix(Ly1) && is.null(Lt1) && is.null(Ymu1) && is.matrix(Ly2) && is.null(Lt2) && is.null(Ymu2)){
    rawCC <- GetRawCrCovFuncFunc(Ly1 = Ly1, Ly2 = Ly2)
    return ( list(smoothedCC = NULL, rawCC = rawCC, bw = NULL, score = NULL) )  
  }
  
  # Otherwise assume you have SPARSE data
  if( is.null(Ymu1) ||   is.null(Ymu2)){
    stop("Both functional means must be provided.")   
  }  
  
  # Get the Raw Cross-covariance    
  rawCC = GetRawCrCovFuncFunc(Ly1 = Ly1, Lt1 = Lt1, Ymu1 = Ymu1, Ly2 = Ly2, Lt2 = Lt2, Ymu2 = Ymu2)
  # Use a heuristic to decide when to bin
  if (!useGAM && sum(duplicated(rawCC$tpairn)) >= 0.2 * length(rawCC$rawCCov)) {
    # message('Binning rawCC')
    rawCC <- BinRawCov(rawCC)
  } else {
    # rename the fields to be the same as the binned rawCC
    names(rawCC)[names(rawCC) == 'rawCCov'] <- 'meanVals'
    names(rawCC)[names(rawCC) == 'tpairn'] <- 'tPairs'
    rawCC$count <- rep(1, length(rawCC$meanVals))
  }
  
  if (rmDiag) {
    diagInd <- rawCC$tPairs[, 1] == rawCC$tPairs[, 2]
    rawCC$tDiag <- rawCC$tPairs[diagInd, , drop=FALSE]
    rawCC$diagMeans <- rawCC$meanVals[diagInd]
    rawCC$diagCount <- rawCC$count[diagInd]
    rawCC$diagRSS <- rawCC$RSS[diagInd]
    rawCC$tPairs <- rawCC$tPairs[!diagInd, , drop=FALSE]
    rawCC$meanVals <- rawCC$meanVals[!diagInd]
    rawCC$count <- rawCC$count[!diagInd]
    rawCC$RSS <- rawCC$RSS[!diagInd]
  }
  
  # Calculate the observation and the working grids
  ulLt1 = unlist(Lt1);             ulLt2 = unlist(Lt2)
  obsGrid1 = sort(unique(ulLt1));  obsGrid2 = sort(unique(ulLt2))
  
  workGrid1 = seq(obsGrid1[1], max(obsGrid1), length.out = 51)
  workGrid2 = seq(obsGrid2[1], max(obsGrid2), length.out = 51)
  workGrid12 = matrix(c(workGrid1, workGrid2),ncol= 2)
  
  if (useGAM == TRUE){ 
    # stop('Cannot be used on binned rawCC') # we cannot use binning for useGAM
    Qdata = data.frame(x =  rawCC$tPairs[,1], y = rawCC$tPairs[,2], z = rawCC$meanVals, group = rawCC$IDs  )
    # I comparsed with 're', ds' and 'gp' too, and 'tp' seems to have a slight edge for what we want
    # myGAM = mgcv::gamm( z ~ s(x,y, bs =c('tp','tp')), random=list(group=~1) , data= Qdata)$gam
    myGAM = mgcv::gam( z ~ s(x,y, bs =c('tp','tp')), data= Qdata)
    estPoints = data.frame( x= rep(workGrid1, times=51), y= rep(workGrid2, each =51), group = rep(3,51*51 )) 
    smoothedCC = matrix(predict(myGAM, estPoints), 51)
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw = NULL, score = NULL) )  
  }
  
  # If the bandwidth is known already smooth the raw CrCov
  if( is.numeric(bw1) &&  is.numeric(bw2)){
    smoothedCC <- smoothRCC2D(rcov =rawCC, bw1, bw2, workGrid1, workGrid2,
                              kern=kern)    
    # potentially incorrect GCV score if the kernel is non-Gaussian
    score = GCVgauss2D(smoothedCC = smoothedCC, smoothGrid = workGrid12, 
                       rawCC = rawCC$meanVals, rawGrid = rawCC$tPairs, 
                       bw1 = bw1, bw2 = bw2)                      
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw =  c(bw1, bw2), score = score, smoothGrid = workGrid12 ) )
    
    # If the bandwidths are unknown use GCV to take find it
  } else {
    # Construct candidate bw's
    bwCandidates <- getBWidths(Lt1, Lt2)
    minGcvScores = Inf
    
    if(bwRoutine == 'grid-search'){
      # Find their associated GCV scores 
      gcvScores = rep(Inf, nrow(bwCandidates)) 
      for (i in 1:length(bwCandidates)){
        gcvScores[i] = theCostFunc(bwCandidates[i,], rawCC, workGrid1, workGrid2, kern, workGrid12)  
      } 
      # Pick the one with the smallest score 
      minimumScore =  min(gcvScores, na.rm=TRUE)
      if( minimumScore == Inf) {
        stop("It seems that the minimum GCV score equals infinity. Stopping 'GetCrCovYX'")
      }
      bInd = which(gcvScores == minimumScore);
      bOpt1 = max(bwCandidates[bInd,1]);
      bOpt2 = max(bwCandidates[bInd,2]); 
      minGcvScores = minimumScore
    } else {
      
      bwRanges = apply(bwCandidates,2, range)
      upperB = bwRanges[2,]
      lowerB = bwRanges[1,]
      
      if( !requireNamespace("minqa", quietly=TRUE) && bwRoutine == 'bobyqa'){ #!is.element('minqa', installed.packages()[,1])
        warning("Cannot use 'minqa::bobyqa' to find the optimal bandwidths. 'minqa' is not installed. We will do an 'L-BFGS-B' search.")
        bwRoutine == 'l-bfgs-b'
      }
      
      if( bwRoutine == 'l-bfgs-b' ){
        theSols = optim(fn = theCostFunc, par = upperB*0.95, # Starting value that "is safe"
                        upper = upperB, lower = lowerB, method ='L-BFGS-B', control = list(maxit = 51),
                        rawCC = rawCC, workGrid1 = workGrid1, workGrid2 = workGrid2, kern = kern, workGrid12 = workGrid12) 
        minGcvScores = theSols$value
      } else { # when BOBYQA is available
        theSols = minqa::bobyqa(fn = theCostFunc, par = upperB*0.95, # Starting value that "is safe"
                                upper = upperB, lower = lowerB, control = list(maxfun = 41),
                                rawCC, workGrid1, workGrid2, kern, workGrid12) 
        minGcvScores = theSols$fval
      }
      
      bOpt1 = theSols$par[1]
      bOpt2 = theSols$par[2]
      
      if( bOpt1 > 0.75 * upperB[1] && bOpt2 > 0.75 * upperB[2] ){
        warning('It seems that the bandwidth selected by the solvers is somewhat large. Maybe you are in local minima.')
      }
    }
    smoothedCC <- smoothRCC2D(rcov=rawCC, bw1 =bOpt1, bw2 =bOpt2, workGrid1, workGrid2, kern=kern)
    return ( list(smoothedCC = smoothedCC, rawCC = rawCC, bw = c(bOpt1, bOpt2), smoothGrid = workGrid12, 
                  score = minGcvScores) )
  }  
}

theCostFunc <- function(xBW, rawCC, workGrid1, workGrid2, kern, workGrid12){
  smoothedCC <- try(silent=TRUE, smoothRCC2D(rcov=rawCC, bw1 = xBW[1], bw2 = xBW[2], 
                                             workGrid1, workGrid2, kern=kern) )
  if( is.numeric(smoothedCC) ){
    theCost = GCVgauss2D( smoothedCC = smoothedCC, smoothGrid = workGrid12, 
                          rawCC = rawCC$meanVals, rawGrid = rawCC$tPairs, bw1 = xBW[1], bw2 = xBW[2])
  } else {
    theCost = Inf
  }
  return(theCost)
}

getBWidths <- function(Lt1, Lt2){
  
  numPoints = 10;
  oldVersion = FALSE
  if(oldVersion == TRUE){
    
 #   bwCandidates <- matrix(rep(0,numPoints * numPoints * 2),ncol=2)
 #   h0 = 2.0 * Minb( sort(ulLt1), 2+1); # 2x the bandwidth needed for at least 3 points in a window
 #   r = diff(range(ulLt1))    
 #   q = (r/(4*h0))^(1/9);  
 #   bwCandidates[,1] = rep( sort(q^( seq(0,12,length.out=numPoints) )*h0), times= numPoints);
 #   h0 = 2.0 * Minb( sort(ulLt2), 2+1); # 2x the bandwidth needed for at least 3 points in a window
 #   r1 = diff(range(ulLt2))    
 #   q = (r/(4*h0))^(1/9);  
 #   bwCandidates[,2] = rep( sort(q^( seq(0,12,length.out=numPoints) )*h0), each= numPoints); 
    
  } else { 
    
    bwCandidates <- matrix(rep(0,numPoints * numPoints * 2),ncol=2) 
    bwCandidates[,1] = rep( seq( BwNN(Lt1, onlyCov = TRUE, k = 5), BwNN(Lt1, onlyCov = TRUE, k = 45), 
                                 length.out=numPoints), times= numPoints); 
    bwCandidates[,2] = rep( seq( BwNN(Lt2, onlyCov = TRUE, k = 5), BwNN(Lt2, onlyCov = TRUE, k = 45), 
                                 length.out=numPoints), each= numPoints); 
    
  } 
  
  return(bwCandidates)
}

smoothRCC2D <- function(rcov,bw1, bw2, xout1, xout2, kern='gauss'){
  # Calculate the smooth Covariances between two functional variables
  # rcov    : raw cross covariance list object returned by BinRawCov(GetRawCrCovFuncFunc())
  # bw1     : scalar
  # bw2     : scalar
  # xout1   : vector M-1
  # xout2   : vector L-1
  # returns : matrix M-L 
  # browser()
  return( Lwls2D( bw = c(bw1, bw2), kern = kern, xin=rcov$tPairs, 
                  yin=rcov$meanVals, win=rcov$count, xout1=xout1, xout2=xout2, crosscov=TRUE) )  
}

GCVgauss2D <- function( smoothedCC, smoothGrid, rawCC, rawGrid, bw1, bw2){ 
  # Calculate GCV cost off smoothed sample assuming a Gaussian kernel
  # smoothedY : vector M-1 
  # smoothedX : vector M-1
  # rawX      : vector N-1
  # rawY      : vector N-1
  # bw        : scalar
  # returns   : scalar
  obsFit <- interp2lin(smoothGrid[,1], smoothGrid[,2], smoothedCC, as.numeric(rawGrid[, 1]), 
                       as.numeric(rawGrid[, 2]))    
  # workaround for degenerate case.
  if (any(is.nan(obsFit))  || any(is.infinite(obsFit))  ){
    return(Inf)
  }  
  # residual sum of squares
  cvsum <- sum((rawCC - obsFit) ^ 2)
  N   = length( rawCC )
  r1  = diff( range(smoothGrid[,1] ) )
  r2  = diff( range(smoothGrid[,2] ) )
  k0 = 0.398942;  # hard-coded constant for Gaussian kernel
  return( cvsum / (1 - (1/N) * (r1 * k0 * r2 * k0) /(bw1 * bw2))^2 )
}
